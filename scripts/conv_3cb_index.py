#!/usr/bin/env python3

"""
This script takes a FASTQ file created by UMI-tools with cell barcodes and UMIs in the header and replaces 
the cell barcodes with an index based on the barcode round files. This index is used to identify the cells
in the downstream analysis, and is identical to the approach BD Rhapsody uses.
The output is two FASTQ files with the new headers.

Usage:
    python conv_3cb_index.py \
        --input_cDNA <input_cDNA FASTQ.gz> \
        --output_cDNA <Output output_cDNA FASTQ.gz> \
        --workdir <Work directory> \
        --bcdir <Barcode directory>

"""

import sys
import os
import gzip
import argparse
import re
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_cDNA", default="/users/asebe/bvanwaardenburg/git/data/241106_BD_Rhapsody_Nvec/Nvec_FINAL/demultiplex/demux_umitools/demux_BCA005_R2.fastq.gz", help="R2 FASTQ.gz (cDNA)")
    parser.add_argument("--output_cDNA", default="/users/asebe/bvanwaardenburg/git/data/241106_BD_Rhapsody_Nvec/Nvec_FINAL/demultiplex/demux_umitools/indexed_BCA005_R2.fastq.gz", help="Output R2 FASTQ.gz")
    parser.add_argument("--workdir", default="/users/asebe/bvanwaardenburg/git/data/241106_BD_Rhapsody_Nvec/Nvec_FINAL/demultiplex/demux_umitools", help="Work directory")
    parser.add_argument("--bcdir", default="/users/asebe/bvanwaardenburg/git/bca_preprocessing/seq_techniques/bd_rhapsody", help="Directory containing barcode files")
    return parser.parse_args()


def split_header_into_CBs(header):
    # Split header into CBs based on underscores
    header = header.strip()
    split_header = header.split("_")
    print(split_header)

    # Assign CBS and UMI
    CB1 = split_header[1][0:9]
    CB2 = split_header[1][9:18]
    CB3 = split_header[1][18:27]
    UMI = split_header[2].split(" ")[0]
    return CB1, CB2, CB3, UMI


def find_CB_in_file(CB, file):
    def hamming_distance(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    # Loops through file to find CB with up to 1 mismatch and returns line
    with open(file, "r") as f:
        for i, line in enumerate(f):
            barcode = line.strip()
            if hamming_distance(CB, barcode) <= 1:
                print("Found CB in line", i)
                return i
        # If it cannot find the CB, return None which will print an error message and skips the read
        else:
            print("CB not found in file")
            return None


def linenum_to_index(CB1_linenum, CB2_linenum, CB3_linenum):
    return (CB1_linenum - 1) * 384 * 384 + (CB2_linenum - 1) * 384 + (CB3_linenum - 1) + 1


def replace_header_with_index(header, index, UMI):
    # Split header into list
    header = header.strip().split(" ")

    # Delete the old cell barcode and UMI assignment
    header[0] = header[0].split("_")[0]

    # Replace the cell barcode with the index, and add the UMI
    header[0] = header[0] + ";" + str(index) + ";" + UMI
    return " ".join(header)


def save_new_fastq(line, o2):
    with gzip.open(o2, "at") as f:
        f.write(line + "\n")


def main():
    args = parse_args()
    os.makedirs(args.bcdir, exist_ok=True)

    # Loop through read 2 and extract the header
    with gzip.open(args.input_cDNA, "rt") as f:
        for line in f:
            if line.startswith("@"):
                header = line
                print(header)

                # Split header into CBs and UMI
                CB1, CB2, CB3, UMI = split_header_into_CBs(header)
                print(CB1, CB2, CB3, UMI)  
                
                # Open barcode round files and search for CB and return line number
                CB1_linenum = find_CB_in_file(CB1, args.bcdir + "/barcodes_R1.txt") 
                CB2_linenum = find_CB_in_file(CB2, args.bcdir + "/barcodes_R2.txt")
                CB3_linenum = find_CB_in_file(CB3, args.bcdir + "/barcodes_R3.txt")
                
                # If one or more CBs are not found, skip the read
                if CB1_linenum is None or CB2_linenum is None or CB3_linenum is None:
                    print("One or more CBs not found in barcode files")
                    continue

                # Converts line number to index
                index = linenum_to_index(CB1_linenum, CB2_linenum, CB3_linenum)
                print("index: ", index)

                # Replace header with index and save to a new FASTQ file
                new_header = replace_header_with_index(header, index, UMI)

                # Write all lines to new FASTQ files
                save_new_fastq(new_header, args.output_cDNA)
            else:
                save_new_fastq(line, args.output_cDNA)

if __name__ == "__main__":
    main()

