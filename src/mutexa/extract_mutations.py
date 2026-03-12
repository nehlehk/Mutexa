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
            if not (str(refline[loc - 1]) == ref):
                print("Error: Ref specified in mutation file doesnt match reference at location: ", loc)
                mismatch_ref = True
        if mismatch_ref:
            print("Error: Try another reference genome!")
            raise SystemExit(1)

    if input_type == 'fasta':
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

    if input_type == 'vcf' or input_type == 'vcfgz':
        positions_to_keep = [int(pos[:-1]) for pos in vcfPosList]
        actual2ref = {i: i for i in range(1, 50000)}
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
            if input_type == 'fasta':
                for i in range(0, len(m), 3):
                    mutPos = m[i]
                    ref = str(m[i + 1])
                    altString = str(m[i + 2])
                    if int(mutPos) in lineDict:
                        snpHaps = categorizeComboHaps(lineDict[int(mutPos)], altString)
                    else:
                        snpHaps = [0] * len(sampleNamelist)
            else:
                for i in range(0, len(m), 3):
                    mutPos, mutRef, mutAlt = m[i:i+3]
                    mutPosInt = int(mutPos)
                    snpHaps = [0] * len(sampleNamelist)
                    if mutPosInt in altDict:
                        vcfAlts = altDict[mutPosInt][0].split(',')
                        alt_index = {vcfAlt: idx + 1 for idx, vcfAlt in enumerate(vcfAlts)}
                        mutAltIndex = alt_index.get(mutAlt, None)
                        if mutAltIndex is not None:
                            lineDictValues = lineDict[mutPosInt][0]
                            for i, item in enumerate(lineDictValues):
                                itemInt = int(item)
                                if itemInt == mutAltIndex:
                                    snpHaps[i] = 1
                                elif itemInt in alt_index.values():
                                    snpHaps[i] = -1
        else:
            try:
                mutPos, mutRef, mutAlt = m
            except:
                print("An exception occurred, the mutation values should be checked")
                continue

            if input_type == 'fasta':
                if int(mutPos) in lineDict:
                    snpHaps = categorizeHap(lineDict[int(mutPos)], mutAlt, mutRef)
                else:
                    snpHaps = [0] * len(sampleNamelist)
            else:
                mutPosInt = int(mutPos)
                snpHaps = [0] * len(sampleNamelist)
                if mutPosInt in altDict:
                    vcfAlts = altDict[mutPosInt][0].split(',')
                    alt_index = {vcfAlt: idx + 1 for idx, vcfAlt in enumerate(vcfAlts)}
                    mutAltIndex = alt_index.get(mutAlt, None)
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




# if __name__ == '__main__':
#     # Define inputs
#     parser = argparse.ArgumentParser(description="Script for extracting variants of interest")
#     parser.add_argument('--mutFile', '-u', help='mutation list (csv)', required=True)
#     parser.add_argument('--alignFile', '-a', help='input alignment file (fasta)/ VCF file', required=True)
#     parser.add_argument('--ref', '-r', help='reference fasta file')
#     parser.add_argument('--prefix','-p', help='outputs perfix')
#     args = parser.parse_args()

#     mutFile = args.mutFile
#     alignFile = args.alignFile
#     refFile = args.ref

#     # Ensure prefix is set correctly
#     if args.prefix is None or args.prefix == '':
#         prefix = "output"
#     else:
#         prefix = args.prefix

#     # mutFile = "ExampleData/mutation.csv"
#     # alignFile = "inputs/VCFs/USA_subset.vcf"
#     # prefix = "VCF"
#     # alignFile = "ExampleData/seqs.fasta"
#     # refFile = "ExampleData/refSeq.fasta"
#     # prefix = "example"

#     # mutFile = "/Users/kar145/0_Work/ClimateChange/data/Dengue_Carol/dengue_mutlist.csv"
#     # alignFile = "/Users/kar145/0_Work/ClimateChange/data/Dengue_Carol/redone_alignment_28022024/dengue_28022024_aligned.fasta"
#     # refFile = "/Users/kar145/0_Work/ClimateChange/data/Dengue_Carol/refKM403575.1_DengueIGIII_WT.fasta"
#     # prefix = "dengue_fasta"
#     # alignFile = "/Users/kar145/0_Work/ClimateChange/data/Dengue_Carol/redone_alignment_28022024/dengue_carol.vcf"
#     # prefix = "dengue_VCF"

#     input_type = check_file_extension(alignFile)

#     #### Carol: Added in for only using the reference if input is fasta, otherwise we can
#     #### read the ref directly from the vcf
#     if input_type == 'fasta':
#         seqlist = list(tqdm(SeqIO.parse(alignFile, 'fasta')))
#     # read reference file
#         if (refFile is None) or (refFile == ""):
#             refline = seqlist[0].seq
#         else:
#             reference = list(SeqIO.parse(refFile, 'fasta'))
#             refline = reference[0].seq
#     else:
#         # vcf_file = pysam.VariantFile(alignFile, "r")
#         print('skipping reference')
#         pass

#     # seqlist = list(tqdm(SeqIO.parse(alignFile, 'fasta')))
#     # # read reference file
#     # if (refFile is None) or (refFile == ""):
#     #     refline = seqlist[0].seq
#     # else:
#     #     reference = list(SeqIO.parse(refFile, 'fasta'))
#     #     refline = reference[0].seq

#     # set outputs
#     mutlist_fname = "outputs/" + prefix + "_mutlist.csv"
#     hapdict_fname = "outputs/" + prefix + "_hapdict.csv"

#     data = pd.read_csv(mutFile, sep='[:]', engine='python', header=None)
#     with open(mutFile, 'r') as fp:
#         for line in enumerate(fp):
#             # search string
#             if ':' in str(line):
#                 data.drop(data.columns[0], axis=1, inplace=True)
#                 break

#     mutList = []
#     for i, row in data.iterrows():
#         line = data.iloc[i, 0]
#         parts = str(line).split(',')
#         mutList.append(parts)

#     # print(mutList)
#     data.to_csv(mutlist_fname, header=False, index=False)


#     vcfPosList = []
#     for m in mutList:
#         for i in range(0, len(m), 3):
#             vcfPosList.append(str(m[i]) + str(m[i + 1]))

#     vcfPosList = sorted(set(vcfPosList))
#     sampleNamelist = []
#     lineDict = defaultdict(list)
#     altDict = defaultdict(list)

#     #### Carol: Similar to above this checks the ref pos and allele in the reference line
#     if input_type == 'fasta':
#         mismatch_ref = False
#         for pos in vcfPosList:
#             ref = pos[-1]
#             loc = int(pos[:-1])
#             if not (str(refline[loc - 1]) == ref):
#                 print("Error: Ref specified in mutation file doesnt match reference at location: ", loc)
#                 mismatch_ref = True
#         if mismatch_ref:
#             print("Error: Try another reference genome!")
#             exit(1)

#     if input_type == 'fasta':
#         with open(alignFile, mode='r') as handle:
#             for record in tqdm(SeqIO.parse(handle, 'fasta')):
#                 sampleNamelist.append(extract_sample_id(record.id))
#                 line = str(record.seq).upper()  # Convert sequence to uppercase
#                 for pos in vcfPosList:
#                     ref = pos[-1].upper()  # Ensure ref is also uppercase
#                     loc = int(pos[:-1])
#                     if str(refline[loc - 1]).upper() == ref:  # Ensure refline comparison is also uppercase
#                         lineDict[loc].append(line[loc - 1])
#                     else:
#                         print("Error: Ref specified in mutation file doesn't match reference at location:", loc)
#                         # exit(1)

#     # print(lineDict)
#     if input_type == 'vcf' or input_type == 'vcfgz':
#         positions_to_keep = [int(pos[:-1]) for pos in vcfPosList]
#         actual2ref = dict()
#         actual2ref = {i: i for i in range(1, 50000)}
#         # with open(alignFile, mode='r') as handle:
#         with (gzip.open if input_type == "vcfgz" else open)(alignFile, 'rt') as handle:
#             for line in handle:
#                 if line.startswith('#CHROM'):
#                     fields = line.strip().split('\t')
#                     sampleNamelist.extend(fields[9:])
#                     continue

#                 if line.startswith('#'):
#                     continue

#                 fields = line.strip().split('\t')
#                 chrom, pos, ref, alt = fields[:4]

#                 if int(pos) in positions_to_keep:
#                     # lineDict[actual2ref[int(pos)]].append(line)
#                     altDict[actual2ref[int(pos)]].append(fields[4])
#                     lineDict[actual2ref[int(pos)]].append(fields[9:])
#     # if input_type == 'vcfgz':
#     #     positions_to_keep = [int(pos[:-1]) for pos in vcfPosList]
#     #     actual2ref = dict()
#     #     actual2ref = {i: i for i in range(1, 50000)}
#     #     with gzip.open(alignFile, mode='rt') as handle: #### Carol: read gzip file
#     #         for line in tqdm(handle):
#     #             # Split the VCF line into fields
#     #             if line.startswith('#CHROM'):
#     #                 fields = line.strip().split('\t')
#     #                 sampleNamelist.extend(fields[9:])
#     #                 continue
#     #
#     #             if line.startswith('#'):
#     #                 continue
#     #             # Extract relevant information from the VCF fields
#     #             fields = line.strip().split('\t')
#     #             chrom, pos, ref, alt = fields[:4]
#     #
#     #             if int(pos) in positions_to_keep:
#     #                 lineDict[actual2ref[int(pos)]].append(line)
#     print("Extracted variants of interest")
#     print("Categorizing variants")
#     hapDict = dict()
#     for m in mutList:
#         # print(m)
#         if len(m) > 3:  # if it's a combination mutation treat it differently. Here 1 = all mutatns present,
#             # 0 = not all mutants present
#             if input_type == 'fasta':
#                 for i in range(0, len(m), 3):
#                     mutPos = m[i]
#                     ref = str(m[i + 1])
#                     altString = str(m[i + 2])
#                     if int(mutPos) in lineDict:
#                         snpHaps = categorizeComboHaps(lineDict[int(mutPos)], altString)
#                     else:
#                         snpHaps = [0] * len(sampleNamelist)
#             else:
#                 ## making hapdict for VCF file -Nehleh
#                 for i in range(0, len(m), 3):
#                     mutPos, mutRef, mutAlt = m[i:i+3]
#                     mutPosInt = int(mutPos)
#                     snpHaps = [0] * len(sampleNamelist)
#                     if mutPosInt in altDict:
#                         vcfAlts = altDict[mutPosInt][0].split(',')
#                         alt_index = {vcfAlt: idx + 1 for idx, vcfAlt in enumerate(vcfAlts)}
#                         mutAltIndex = alt_index.get(mutAlt, None)
#                         if mutAltIndex is not None:
#                             lineDictValues = lineDict[mutPosInt][0]
#                             for i, item in enumerate(lineDictValues):
#                                 itemInt = int(item)
#                                 if itemInt == mutAltIndex:
#                                     snpHaps[i] = 1
#                                 elif itemInt in alt_index.values():
#                                     snpHaps[i] = -1
#         else:  # if it's just a single SNP: 1 = mutant present, 0 = reference present, -1 = something else present (
#             # other haploytpe or ambiguous base)
#             try:
#                 mutPos, mutRef, mutAlt = m
#             except:
#                 print("An exception occurred, the mutation valuse should be checked")
#             if input_type == 'fasta':
#                 if int(mutPos) in lineDict:  # check if variant was actually found at location in vcf file
#                     snpHaps = categorizeHap(lineDict[int(mutPos)], mutAlt, mutRef)
#                 else:  # if not, set all haplotypes to reference (0)
#                     snpHaps = [0] * len(sampleNamelist)
#             else:
#                 ## making hapdict for VCF file -Nehleh
#                 mutPosInt = int(mutPos)
#                 snpHaps = [0] * len(sampleNamelist)
#                 if mutPosInt in altDict:
#                     vcfAlts = altDict[mutPosInt][0].split(',')
#                     alt_index = {vcfAlt: idx + 1 for idx, vcfAlt in enumerate(vcfAlts)}
#                     mutAltIndex = alt_index.get(mutAlt, None)

#                     if mutAltIndex is not None:
#                         lineDictValues = lineDict[mutPosInt][0]

#                         for i, item in enumerate(lineDictValues):
#                             itemInt = int(item)
#                             if itemInt == mutAltIndex:
#                                 snpHaps[i] = 1
#                             elif itemInt in alt_index.values():
#                                 snpHaps[i] = -1


#         for i in range(0, len(sampleNamelist)):
#             if sampleNamelist[i] in hapDict:
#                 hapDict[sampleNamelist[i]].append(snpHaps[i])
#             else:
#                 hapDict[sampleNamelist[i]] = [snpHaps[i]]

#     df = pd.DataFrame.from_dict(hapDict)
#     final = df.T.reset_index()
#     final.to_csv(hapdict_fname, header=False, index=False)
#     print("Done")

#     # vcf_file.close()

