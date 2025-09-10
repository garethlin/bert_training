import os
import sys
import csv
import configparser
from pathlib import Path
import pandas as pd
#from neo4j import GraphDatabase
import os
import json
from itertools import combinations 
import random
import numpy

def ecnum2class(ecnum):
    return ".".join(ecnum.split(".")[:3])

def connect(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    driver = GraphDatabase.driver(
        config["neo4j"]["host"],
        auth=(
            config["neo4j"]["username"],
            config["neo4j"]["password"],
        ),
    )

    driver.verify_connectivity()
    driver.verify_authentication()

    return driver

def run_query(driver, query, params={}):
    with driver.session() as session:
        result = session.run(query, params)
        return pd.DataFrame([r.values() for r in result], columns=result.keys())

def get_enzymes(driver):
    QUERY = "MATCH (e: Enzyme) RETURN elementID(e) as id, e.sequence as sequence"
    enzymes = run_query(driver, QUERY)
    return enzymes

def get_fasta(enzymes, outfilename):
    with open(outfilename, "w") as fasta_file:
        for _, row in enzymes.iterrows():
            if row['sequence'] != '':
                fasta_file.write(f">{row['id']}\n")
                fasta_file.write(f"{row['sequence']}\n")


def filt_blast_p(query_prots, target_prots, min_len, min_id):
    """
    Perform BLAST search and filter results based on length and identity.
    """
    os.system(f"makeblastdb -in {target_prots} -out {target_prots} -dbtype prot")
    outfmt = "6 qseqid qlen qstart qend sstrand sseqid slen sstart send length pident mismatch gapopen evalue bitscore"
    stdout = os.popen(f'blastp -db {target_prots} -query {query_prots} -outfmt "{outfmt}"').read().rstrip()
    filtered_hits = []
    for raw_hit in stdout.splitlines():
        fields = raw_hit.split()
        qlen = int(fields[1])
        align_len = int(fields[3]) - int(fields[2])
        pid = float(fields[10])
        if align_len / qlen > min_len and pid > min_id:
            filtered_hits.append(raw_hit)
    return filtered_hits


def split_ec(input_file, output_file):
    """
    Process a TSV file with uniprot accession, EC number, and sequence fields.
    Expand records with multiple EC numbers into multiple individual records.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Get the header
        header = next(reader)
        writer.writerow(header)

        # Process each row
        for row in reader:
            if len(row) != 3:
                print(f"Warning: Skipping invalid row (expected 3 columns, found {len(row)}): {row}")
                continue

            uniprot_acc, ec_numbers, sequence = row

            # Split EC numbers by semicolon and strip whitespace
            ec_list = [ec.strip() for ec in ec_numbers.split(';')]

            # Create a new row for each EC number
            for ec in ec_list:
                if ec:  # Skip empty EC numbers
                    writer.writerow([uniprot_acc, ec, sequence])

def tab_to_fasta(input_file, output_file):
    # Read the input file
    with open(input_file, 'r') as file:
        lines = file.read().strip().split('\n')
    
    fasta_content = []
    seq_lookup = {}
    
    for line in lines:
        if not line.strip():
            continue
        if "Sequence" in line:
            continue
            
        parts = line.split('\t')
        if len(parts) < 3:
            continue
            
        accession_id, sequence = parts[0], parts[2]
        seq_lookup[accession_id] = sequence
        
        # Add the FASTA header line
        fasta_content.append(f">{accession_id}")
        
        # Add the entire sequence on a single line
        fasta_content.append(sequence.strip())
    
    fasta_result = '\n'.join(fasta_content)
    
    # Write to file or return the string
    with open(output_file, 'w') as file:
        file.write(fasta_result)
    return seq_lookup


def get_ec_data(input_file):
    data = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            if "EC number" in line:
                next
            else:
                acc = line.rstrip().split()[0]
                ec_number = line.rstrip().split()[1]
                ec_class = ecnum2class(ec_number)
                if "-" in ec_number:
                    pass
                else:
                    if ec_class not in data.keys():
                        data[ec_class] = {}
                    if ec_number not in data[ec_class].keys():
                        data[ec_class][ec_number] = []
                    data[ec_class][ec_number].append(acc)

    non_unit_ec = {}
    by_count = {}

    for ec_class in data:
        for ec_number in data[ec_class]:
            if len(data[ec_class][ec_number]) == 1:
                pass
            else:
                if ec_class not in non_unit_ec.keys():
                    non_unit_ec[ec_class] = {}
                if ec_number not in non_unit_ec[ec_class].keys():
                    non_unit_ec[ec_class][ec_number] = data[ec_class][ec_number]

    return non_unit_ec

def get_accession_lookup(data):
    lookup = {}
    for ec_class in data.keys():
        for ec_number in data[ec_class].keys():
            for acc in data[ec_class][ec_number]:
                lookup[acc] = ec_number
    return lookup
            

def blast2exclude(blast_hits):
    exclude = []
    for hit in blast_hits:
        acc = hit.split()[5]
        exclude.append(acc)
    return exclude

def get_long_protein_headers(fasta_path, min_length=2048):
    """
    Parses a FASTA file and returns a list of headers for proteins longer than min_length amino acids.
    The header returned is the substring up to the first whitespace after the '>'.

    Args:
        fasta_path (str): Path to the FASTA file.
        min_length (int, optional): Minimum protein length to filter. Defaults to 2048.

    Returns:
        List[str]: List of FASTA headers (without '>') for qualifying proteins, up to the first whitespace.
    """
    headers = []
    current_header = None
    current_seq = []

    def flush():
        if current_header is not None:
            seq = ''.join(current_seq)
            if len(seq) > min_length:
                headers.append(current_header)

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                flush()
                # Take the part up to the first whitespace
                current_header = line[1:].strip().split(None, 1)[0]
                current_seq = []
            else:
                current_seq.append(line)
        flush()  # Handle last entry

    return headers

def filter_tsv(exclusion_list, input_file_path, output_file_path):
    """
    Filter a TSV file to exclude rows where the Entry field matches any element in the exclusion list.
    
    Args:
        exclusion_list (list): List of Entry codes to exclude
        input_file_path (str): Path to the input TSV file
        output_file_path (str): Path to write the filtered output TSV
    
    Returns:
        str: Path to the output file
    """
    # Create a set from the exclusion list for faster lookups
    exclusion_set = set(exclusion_list)
    
    try:
        with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
            # Read and write the header
            header = infile.readline()
            outfile.write(header)
            
            # Process each line
            for line in infile:
                fields = line.strip().split('\t')
                
                # Skip if there aren't enough fields
                if len(fields) < 3:
                    continue
                
                # The Entry is the first field (index 0)
                entry = fields[0]
                
                # Check if the Entry field is in the exclusion list
                if entry not in exclusion_set:
                    outfile.write(line)
        
        return output_file_path
    
    except Exception as e:
        raise Exception(f"Error filtering TSV file: {e}")

def select_close_pairs_from_acc_list(accession_list):
    all_pairs = list(combinations(accession_list,2))
    if len(accession_list) > 5:
        selected_pairs = random.sample(all_pairs, 10)
    else:
        selected_pairs = all_pairs
    
    randomized_pairs = []
    for a, b in selected_pairs:
        if random.choice([True, False]):
            randomized_pairs.append((b, a))
        else:
            randomized_pairs.append((a, b))

    assert len(randomized_pairs) != 2

    if len(randomized_pairs) == 1:
        outcomes = ['train', 'test', 'val']
        probabilities = [0.8, 0.1, 0.1]
        choice = numpy.random.choice(outcomes, p=probabilities)
        if choice == 'train':
            train = randomized_pairs
            test = None
            val = None
        elif choice == 'test':
            test = randomized_pairs
            train = None
            val = None
        elif choice == 'val':
            val = randomized_pairs
            test = None
            train = None
    else:
            val = []
            test = []
            train = []
            for pair in randomized_pairs:
                outcomes = ['train', 'test', 'val']
                probabilities = [0.8, 0.1, 0.1]
                choice = numpy.random.choice(outcomes, p=probabilities)
                if choice == 'train':
                    train.append(pair)
                elif choice == 'test':
                    test.append(pair)
                elif choice == 'val':
                    val.append(pair)

    return test, val, train
        
def create_far_pair(pair, data, acc_lookup):
    acc_number = pair[0]
    ec_number = acc_lookup[acc_number]
    ec_class = ecnum2class(ec_number)
    ec_classes = list(data.keys())
    ec_classes.remove(ec_class)
    random_class = random.choice(ec_classes)
    ec_numbers = list(data[random_class].keys())
    random_ec_number = random.choice(ec_numbers)
    accessions = data[random_class][random_ec_number]
    random_accession1 = random.choice(data[ec_class][ec_number])
    random_accession2 = random.choice(data[random_class][random_ec_number])
    return (random_accession1,random_accession2)
    

def create_near_pair(pair):
    ec_number = acc_lookup[pair[0]]
    ec_class = ecnum2class(ec_number)
    ec_number_list = list(data[ec_class].keys())
    if len(ec_number_list) > 1:
        ec_number_list.remove(ec_number)
        random_ec_number = random.choice(ec_number_list)
        random_accession1 = random.choice(data[ec_class][ec_number])
        random_accession2 = random.choice(data[ec_class][random_ec_number])
        return (random_accession1,random_accession2)
    else:
        return None

def get_same_pairs(data):
    test_pairs = []
    val_pairs = []
    train_pairs = []
    for ec_class in data.keys():
        for ec_number in data[ec_class].keys():
            test, val, train = select_close_pairs_from_acc_list(data[ec_class][ec_number])
            if test is not None:
                test_pairs.extend(test)
            if val is not None:
                val_pairs.extend(val)
            if train is not None:
                train_pairs.extend(train)
    return test_pairs, val_pairs, train_pairs

def get_pairs(data, acc_lookup):
    test_pairs, val_pairs, train_pairs = get_same_pairs(data)
    diff_ec_pairs = {}
    diff_ec_pairs['train'] = []
    diff_ec_pairs['test'] = []
    diff_ec_pairs['val'] = []

    same_ec_pairs = {}
    same_ec_pairs['train'] = train_pairs
    same_ec_pairs['test'] = test_pairs
    same_ec_pairs['val'] = val_pairs
    
    for idx, pair_list in enumerate([test_pairs, val_pairs, train_pairs]):
        if idx == 0:
            loop = 'test'
        elif idx == 1:
            loop = 'val'
        elif idx == 2:
            loop = 'train'
        for pair in pair_list:
            far_pair = create_far_pair(pair, data, acc_lookup)
            near_pair = create_near_pair(pair)
            diff_ec_pairs[loop].append(far_pair)
            if near_pair is not None:
                diff_ec_pairs[loop].append(near_pair)
    return same_ec_pairs, diff_ec_pairs

db_config = sys.argv[1]
outdir = sys.argv[2]

os.system(f"mkdir -p {outdir} && curl \"https://rest.uniprot.org/uniprotkb/stream?download=true&fields=accession%2Cec%2Csequence&format=tsv&query=%28%28reviewed%3Atrue%29+AND+%28ec%3A*%29+AND+%28taxonomy_id%3A33090%29%29\" > {outdir}/swissprot_plant_with_ec.tsv")
driver = connect(sys.argv[1])
enzymes = get_enzymes(driver)
enz_fasta_file = f"{outdir}/enzymes.fasta"
sp_fasta_file = f"{outdir}/swissprot.fasta"
get_fasta(enzymes, enz_fasta_file)
seq_lookup = tab_to_fasta(f"{outdir}/swissprot_plant_with_ec.tsv", sp_fasta_file)
blast_hits = filt_blast_p(enz_fasta_file, sp_fasta_file, 0.9, 97)
exclude_list = blast2exclude(blast_hits)
long_prot_list = get_long_protein_headers(sp_fasta_file)
exclude_list.extend(long_prot_list)
filter_tsv(exclude_list, f"{outdir}/swissprot_plant_with_ec.tsv", f"{outdir}/swissprot_plant_with_ec.exclude_graph.tsv")
split_ec(f"{outdir}/swissprot_plant_with_ec.exclude_graph.tsv",f"{outdir}/swissprot_plant_with_ec.exclude_graph.split.tsv")
data = get_ec_data(f"{outdir}/swissprot_plant_with_ec.exclude_graph.split.tsv")
acc_lookup = get_accession_lookup(data)
same_ec_pairs, diff_ec_pairs = get_pairs(data, acc_lookup)

for dset in ['test', 'val', 'train']:
    with open(f"{dset}.pairs.csv", 'w') as ec_check_file:       
        with open(f"{dset}.model_data.csv", 'w') as model_data_file:
            for pair in same_ec_pairs[dset]:
               ec_check_file.write(",".join([pair[0], pair[1], acc_lookup[pair[0]], acc_lookup[pair[1]]]) + "\n")
               model_data_file.write(",".join([ seq_lookup[pair[0]], seq_lookup[pair[1]], "1" ]) + "\n")
            for pair in diff_ec_pairs[dset]:
                ec_check_file.write(",".join([pair[0], pair[1], acc_lookup[pair[0]], acc_lookup[pair[1]]]) + "\n") 
                model_data_file.write(",".join([ seq_lookup[pair[0]], seq_lookup[pair[1]], "0" ]) + "\n")

