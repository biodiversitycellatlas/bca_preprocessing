#!/usr/bin/env python3

import os
import gzip
import argparse
from collections import defaultdict


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return "".join(comp.get(base, base) for base in reversed(seq))


def hamming(s1, s2):
    """
    Compute the Hamming distance between two equal-length strings.
    If lengths differ, returns a large number to avoid a false match.
    """
    if len(s1) != len(s2):
        return float('inf')
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def process_fastq_for_barcodes(fastq_path, barcodes):
    """
    Process a gzipped FASTQ file to assign each read a barcode label.
    
    For each record (assuming standard FASTQ format), the sequence is compared
    to the list of provided barcodes (allowing up to 1 mismatch). The first barcode
    that meets the criterion is assigned, otherwise the read is labeled as 'unfiltered'.
    
    Returns:
        A dictionary mapping read_id -> matched barcode label.
    """
    mapping = {}
    with gzip.open(fastq_path, "rt") as f:
        line_count = 0
        for line in f:
            
            if line_count % 4 == 0:  # header line
                header = line.rstrip("\n")
                read_id = header[1:].split()[0] if header.startswith("@") else header.split()[0]
            
            elif line_count % 4 == 1:  # sequence line
                seq = line.rstrip("\n")
                matched = None
                for bc in barcodes:
                    if hamming(seq, bc) <= 1:
                        matched = bc
                        break
                mapping[read_id] = matched if matched is not None else "unfiltered"
            
            line_count += 1
    return mapping


def extract_fastq_subset(input_fastq, combined_mapping, read_label, output_dir):
    """
    Extract FASTQ records from a gzipped input file for reads found in combined_mapping.
    
    For each read in the input FASTQ file, if its read_id exists in combined_mapping,
    the record is written to an output file named "{barcodeCombo}_{read_label}.fastq.gz"
    in the given output directory.
    """
    # Create file handles for each unique barcode combination.
    output_handles = {}
    for combo in set(combined_mapping.values()):
        out_filename = os.path.join(output_dir, f"{combo}_{read_label}.fastq.gz")
        output_handles[combo] = gzip.open(out_filename, "wt")
    
    with gzip.open(input_fastq, "rt") as f:
        record = []
        for line in f:
            record.append(line)
            if len(record) == 4:
                header = record[0].rstrip("\n")
                read_id = header[1:].split()[0] if header.startswith("@") else header.split()[0]
                if read_id in combined_mapping:
                    combo = combined_mapping[read_id]
                    output_handles[combo].writelines(record)
                record = []
    
    for handle in output_handles.values():
        handle.close()


def main():
    parser = argparse.ArgumentParser(
        description="Split FASTQ files based on i5 and i7 barcode matching (with up to 1 mismatch) and using the reverse-complement of provided barcodes."
    )
    parser.add_argument('--work-dir', type=str,
                        default="/users/asebe/bvanwaardenburg/git/data/250331_OAKseq_Nvec",
                        help="Base working directory.")
    parser.add_argument('--fastq-dir', type=str,
                        default=None,
                        help="Directory containing FASTQ files. Defaults to work_dir/fastq.")
    parser.add_argument('--output-dir', type=str,
                        default=None,
                        help="Directory to store output FASTQs. Defaults to work_dir/fastq_splitted.")
    parser.add_argument('--i5-barcodes', type=str,
                        default="CTGTCCTGCT,AGCGGGATTT,CTTGATCGTA,CGTACCGTTA,GTAGGAGTCG",
                        help="Comma-separated list of i5 barcode sequences (provided in forward orientation).")
    parser.add_argument('--i7-barcodes', type=str,
                        default="TTTGATCCAC,TTTGATCCAC,TTTGATCCAC,TTTGATCCAC,TGATGATTCA",
                        help="Comma-separated list of i7 barcode sequences (provided in forward orientation).")
    parser.add_argument('--R1', type=str,
                        default="Undetermined_S0_R1_001.fastq.gz",
                        help="Filename for R1 FASTQ file.")
    parser.add_argument('--R2', type=str,
                        default="Undetermined_S0_R2_001.fastq.gz",
                        help="Filename for R2 FASTQ file.")
    parser.add_argument('--I2', type=str,
                        default="Undetermined_S0_I2_001.fastq.gz",
                        help="Filename for I2 FASTQ file (i5 indices).")
    parser.add_argument('--I1', type=str,
                        default="Undetermined_S0_I1_001.fastq.gz",
                        help="Filename for I1 FASTQ file (i7 indices).")
    args = parser.parse_args()

    # Set directories
    work_dir = args.work_dir
    fastq_dir = args.fastq_dir if args.fastq_dir else os.path.join(work_dir, "fastq")
    output_dir = args.output_dir if args.output_dir else os.path.join(work_dir, "fastq_splitted")
    os.makedirs(output_dir, exist_ok=True)

    # Parse and reverse-complement barcodes
    i5_barcodes = [reverse_complement(bc.strip()) for bc in args.i5_barcodes.split(",")]
    i7_barcodes = [bc.strip() for bc in args.i7_barcodes.split(",")]

    print("Using reverse-complemented i5 barcodes:", i5_barcodes)
    print("Using i7 barcodes:", i7_barcodes)

    # Define FASTQ file paths
    R1_file = os.path.join(fastq_dir, args.R1)
    R2_file = os.path.join(fastq_dir, args.R2)
    I2_file = os.path.join(fastq_dir, args.I2)
    I1_file = os.path.join(fastq_dir, args.I1)

    # Process barcode assignments
    print("Processing I2 for i5 barcodes...")
    i5_mapping = process_fastq_for_barcodes(I2_file, i5_barcodes)
    print("Processing I1 for i7 barcodes...")
    i7_mapping = process_fastq_for_barcodes(I1_file, i7_barcodes)

    # Build a combined mapping for reads present in both I2 and I1
    combined_mapping = {}
    for read_id in set(i5_mapping) & set(i7_mapping):
        combined_mapping[read_id] = f"{i5_mapping[read_id]}_{i7_mapping[read_id]}"

    # Report counts by barcode combination
    combination_groups = defaultdict(int)
    for combo in combined_mapping.values():
        combination_groups[combo] += 1
    for combo, count in combination_groups.items():
        print(f"Found {count} reads for combination {combo}")

    # Extract corresponding FASTQ subsets
    for read_label, fastq_file in [("R1", R1_file), ("R2", R2_file), ("I2", I2_file), ("I1", I1_file)]:
        print(f"Processing {read_label} file: {fastq_file}")
        extract_fastq_subset(fastq_file, combined_mapping, read_label, output_dir)

    print("Done. Split FASTQs are in:", output_dir)

if __name__ == "__main__":
    main()
