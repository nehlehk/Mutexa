import warnings
import csv
import pandas as pd
import re
import argparse




def check_format(line):
    # Define a case-insensitive regular expression pattern for the expected format
    mutation_pattern = re.compile(r'(?i)[-A-Z0-9_.+]+\:\d+,[ACGT\*],[-ACGT\*](,\d+,[ACGT\*],[-ACGT\*])*')

    # Check if the line contains at least one occurrence of the mutation pattern
    if mutation_pattern.fullmatch(line):
        return True
    else:
        return False

def check_mutfile_format(mutFile):
    mutation_ids_set = set()
    with open(mutFile, 'r') as file:
        for line in file:
            line = line.strip()
            if not check_format(line):
                print(f"Error: Incorrect format in line - {line}")
                exit(1)

            mutation_id = line.split(":")[0].lower()  # Extract and convert to lowercase for case-insensitivity

            if mutation_id in mutation_ids_set:
                print(f"Error: Duplicate mutation_id found - {mutation_id}")
                exit(1)

            mutation_ids_set.add(mutation_id)
def check_metafile_format(metaFile):
    df = pd.read_csv(metaFile)
    print("Checking the required columns are available ... ")
    flag_id = False
    flag_date = False
    flag_loc = False
    for name in df.columns:
        if re.search(r'accession|id|SampleID', name, re.IGNORECASE):
            flag_id = True
        if re.search(r'collection|date|date_submitted|SampleCollectionDate|date', name, re.IGNORECASE):
            flag_date = True
            date_str = str(name)
        if re.search(r'continent|country|region|location|division|lineage|group|strain', name, re.IGNORECASE):
            flag_loc = True

    if not flag_id:
        print("Error: The CSV file is missing a column related to ID.")
        exit(1)
    if not flag_date:
        print("Error: The CSV file is missing a column related to collection date.")
        exit(1)
    if not flag_loc:
        print("Error: The CSV file is missing a column related to location.")
        exit(1)


warnings.filterwarnings("ignore")
if __name__ == '__main__':
    # Define inputs
    parser = argparse.ArgumentParser(description='Merge metadata with mutation ')
    parser.add_argument('--mutFile', '-u', help='mutation list (csv)', required=True)
    parser.add_argument('--metaFile','-m', help = 'Input metadata file (.json or .csv)',required=True)
    # parser.add_argument('--ref', '-r', help='reference fasta file')


    args = parser.parse_args()
    mutFile = args.mutFile
    metaFile = args.metaFile
    # refFile = args.ref

    # metaFile = '/Users/kar145/Projects/strepifun/inputs/testing_files/MPXV-test-data/metadata031123_subset-all.csv'
    # mutFile = '/Users/kar145/Projects/strepifun/inputs/testing_files/MPXV-test-data/mpxv_mutlist_new.csv'

    # check the format of mutation file
    check_mutfile_format(mutFile)
    # check the format of metadata file
    check_metafile_format(metaFile)

