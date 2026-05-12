import time
from Bio import SeqIO
from tqdm import tqdm
import argparse
from collections import defaultdict
import csv
import pandas as pd
import re
import os
import gzip
from typing import Optional
# import pysam

# Creates mutlist and hapdict file for use in merge_meta.py

def check_file_extension(file_path):
    # Get the file extension
    _, file_extension = os.path.splitext(file_path)
    # Check if the file extension is either 'fasta' or 'vcf'
    if (file_extension.lower() == '.fasta') or (file_extension.lower() == '.fa'):
        return 'fasta'
    elif file_extension.lower() == '.vcf':
        return 'vcf'
    elif file_extension.lower() == '.gz':
        return 'vcfgz'
    else:
        return None


def extract_sample_id(string):
    match = re.search(r'\|([^|]+)\|', string)
    if match:
        sample_id = match.group(1)
        return sample_id
    else:
        match = re.search(r'([^/|]+)', string)
        if match:
            sample_id = match.group(1)
            return sample_id
def categorizeHap(altList, altString, ref):
    # this function categorizes the haplotypes of each sample as ref (0), alt (1) or other (-1)
    cat = []
    for hap in altList:
        if hap == ref:
            cat.append(0)
        elif hap == altString:
            cat.append(1)
        else:
            cat.append(-1)
    return cat


def categorizeComboHaps(altList, altString):
    # haps, alts
    cat = []
    for i in range(0, len(altList)):
        if altList[i] == altString:
            cat.append(1)
        else:
            cat.append(0)
    return cat


def run_extract(mut_file: str, align_file: str, ref_genome: Optional[str] = None, prefix: str = "output") -> None:
    """
    Extract variants of interest and write:
      outputs/<prefix>_mutlist.csv
      outputs/<prefix>_hapdict.csv
    """

    # ---- BEGIN: moved from your __main__ block ----
    mutFile = mut_file
    alignFile = align_file
    refFile = ref_genome

    # Ensure prefix is set correctly
    if prefix is None or prefix == "":
        prefix = "output"

    input_type = check_file_extension(alignFile)

    if input_type == 'fasta':
        seqlist = list(tqdm(SeqIO.parse(alignFile, 'fasta')))
        # read reference file
        if (refFile is None) or (refFile == ""):
            refline = seqlist[0].seq
        else:
            reference = list(SeqIO.parse(refFile, 'fasta'))
            refline = reference[0].seq
    else:
        print('skipping reference')
        # (VCF case) refline unused below except for fasta checks

    # set outputs
    os.makedirs("outputs", exist_ok=True)
    mutlist_fname = "outputs/" + prefix + "_mutlist.csv"
    hapdict_fname = "outputs/" + prefix + "_hapdict.csv"

    data = pd.read_csv(mutFile, sep='[:]', engine='python', header=None)
    with open(mutFile, 'r') as fp:
        for line in enumerate(fp):
            if ':' in str(line):
                data.drop(data.columns[0], axis=1, inplace=True)
                break

    mutList = []
    for i, row in data.iterrows():
        line = data.iloc[i, 0]
        parts = str(line).split(',')
        mutList.append(parts)

    data.to_csv(mutlist_fname, header=False, index=False)

    vcfPosList = []
    for m in mutList:
        for i in range(0, len(m), 3):
            vcfPosList.append(str(m[i]) + str(m[i + 1]))

    vcfPosList = sorted(set(vcfPosList))
    sampleNamelist = []
    lineDict = defaultdict(list)
    altDict = defaultdict(list)

    if input_type == 'fasta':
        mismatch_ref = False
        for pos in vcfPosList:
            ref = pos[-1]
            loc = int(pos[:-1])
            if not (str(refline[loc - 1]).upper() == ref):
                print("Error: Ref specified in mutation file doesnt match reference at location: ", loc)
                mismatch_ref = True
        if mismatch_ref:
            print("Error: Try another reference genome!")
            raise SystemExit(1)

    if input_type == 'fasta':
        print("Reading FASTA file...")

        with open(alignFile, mode='r') as handle:
            for record in tqdm(SeqIO.parse(handle, 'fasta')):
                sampleNamelist.append(extract_sample_id(record.id))
                line = str(record.seq).upper()
                for pos in vcfPosList:
                    ref = pos[-1].upper()
                    loc = int(pos[:-1])
                    if str(refline[loc - 1]).upper() == ref:
                        lineDict[loc].append(line[loc - 1])
                    else:
                        print("Error: Ref specified in mutation file doesn't match reference at location:", loc)
    sampleNamelist = list(set(sampleNamelist))
    if input_type == 'vcf' or input_type == 'vcfgz':
        print("Reading VCF...")
        positions_to_keep = [int(pos[:-1]) for pos in vcfPosList]
        actual2ref = {i: i for i in range(1, 50000)} # max genome of 50k, increase if genome is larger
        with (gzip.open if input_type == "vcfgz" else open)(alignFile, 'rt') as handle:
            for line in handle:
                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    sampleNamelist.extend(fields[9:])
                    continue
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                chrom, pos, ref, alt = fields[:4]
                if int(pos) in positions_to_keep:
                    altDict[actual2ref[int(pos)]].append(fields[4])
                    lineDict[actual2ref[int(pos)]].append(fields[9:])

    print("Extracted variants of interest")
    print("Categorizing variants")

    hapDict = dict()
    for m in mutList:
        if len(m) > 3:
            sample_altHaps = [] # sample mutation set
            altHaps = [] # mutation of interest index in VCF/FASTA set
            mutation_pos = [] # positions to subset alt library
            
            if input_type == 'fasta':
                m = [item.replace("*", "-") for item in m] # Replace mutlist deletion character if vcf coded
                
                for i in range(0, len(m), 3):
                    mutPos = m[i]
                    ref = str(m[i + 1])
                    altString = str(m[i + 2])
                    mutPosInt = int(mutPos) 
                    mutation_pos.append(mutPosInt) # Add to list of mutation positions for mutation set of interest
                    altHaps.append(altString) # Add set of mutation set of interest

                if altHaps is not None:
                    temp_dict = {pos: lineDict[pos] for pos in mutation_pos if pos in lineDict} # Subset mutation set
                    sample_altHaps = list(zip(*temp_dict.values()))
                    sample_altHaps= [list(haps) for haps in sample_altHaps] 

                    snpHaps = categorizeComboHaps(sample_altHaps, altHaps)
                else:
                    snpHaps = [0] * len(sampleNamelist)

            else:
                m = [item.replace("-", "*") for item in m] # Replace mutlist deletion character if fasta coded

                for i in range(0, len(m), 3):
                    mutPos, mutRef, mutAlt = m[i:i+3]
                    mutPosInt = int(mutPos)
                    snpHaps = [0] * len(sampleNamelist)
                    mutation_pos.append(mutPosInt)

                    if mutPosInt in altDict:
                        vcfAlts = altDict[mutPosInt][0].split(',')
                        alt_index = {vcfAlt: idx + 1 for idx, vcfAlt in enumerate(vcfAlts)}
                        altHaps.append(alt_index.get(mutAlt))
                        if altHaps is not None:
                            temp_dict = {pos: lineDict[pos] for pos in mutation_pos if pos in lineDict}
                            flattened_temp_dict = {k: [item for sublist in v for item in sublist] for k, v in temp_dict.items()}
                            
                            for i in range(0,int(len(sampleNamelist))):
                                sample_altHaps=[int(sub[i]) for sub in flattened_temp_dict.values()]
                                if sample_altHaps == altHaps:
                                    snpHaps[i] = 1
                                    #print(sample_altHaps,altHaps)
                                    #snpHaps[i]=(sample_altHaps,altHaps) # FOr checking sample alts and althaps
                                elif sample_altHaps!= altHaps:
                                    snpHaps[i] = 0
                        
        else:
            try:
                mutPos, mutRef, mutAlt = m
            except:
                print("An exception occurred, the mutation values should be checked")
                continue

            if input_type == 'fasta':
                m = [item.replace("*", "-") for item in m] # Replace mutlist deletion character if vcf coded
                mutPos, mutRef, mutAlt = m
                if int(mutPos) in lineDict:
                    snpHaps = categorizeHap(lineDict[int(mutPos)], mutAlt, mutRef)
                else:
                    snpHaps = [0] * len(sampleNamelist)
            else:
                m = [item.replace("-", "*") for item in m]  # Replace mutlist deletion character if fasta coded
                mutPos, mutRef, mutAlt = m
                mutPosInt = int(mutPos)
                snpHaps = [0] * len(sampleNamelist)
                if mutPosInt in altDict:
                    vcfAlts = altDict[mutPosInt][0].split(',')
                    alt_index = {vcfAlt: idx + 1 for idx, vcfAlt in enumerate(vcfAlts)}
                    mutAltIndex = alt_index.get(mutAlt, -1)
                    if mutAltIndex is not None:
                        lineDictValues = lineDict[mutPosInt][0]
                        for i, item in enumerate(lineDictValues):
                            itemInt = int(item)
                            if itemInt == mutAltIndex:
                                snpHaps[i] = 1
                            elif itemInt in alt_index.values():
                                snpHaps[i] = -1

        for i in range(0, len(sampleNamelist)):
            if sampleNamelist[i] in hapDict:
                hapDict[sampleNamelist[i]].append(snpHaps[i])
            else:
                hapDict[sampleNamelist[i]] = [snpHaps[i]]

    df = pd.DataFrame.from_dict(hapDict)
    final = df.T.reset_index()
    final.to_csv(hapdict_fname, header=False, index=False)
    print("Done")
    # ---- END ----


def main():
    parser = argparse.ArgumentParser(description="Script for extracting variants of interest")
    parser.add_argument('--mutFile', '-u', help='mutation list (csv)', required=True)
    parser.add_argument('--alignFile', '-a', help='input alignment file (fasta)/ VCF file', required=True)
    parser.add_argument('--ref', '-r', help='reference fasta file')
    parser.add_argument('--prefix','-p', help='outputs prefix')
    args = parser.parse_args()

    run_extract(
        mut_file=args.mutFile,
        align_file=args.alignFile,
        ref_genome=args.ref,
        prefix=args.prefix or "output",
    )


if __name__ == '__main__':
    main()