#!/usr/bin/env python
"""
A fastq splitter for single cell RNA-seq data.

This script:
  - Loads a barcode whitelist from a JSON file (mapping barcode -> group)
  - Reads paired fastq (R1 and R2) records from gzip-compressed files
  - Extracts a barcode from a fixed position in the R2 sequence
  - Checks if the barcode exists in the whitelist and, if so,
    writes the record pair to group-specific output files.
"""

import argparse
import gzip
import os

def load_whitelist_csv(file_path):
    """Load a CSV whitelist mapping barcode sequence to well."""
    whitelist = {}
    with open(file_path, 'r') as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 5:
                continue
            bci, sequence, uid, well, stype = parts
            whitelist[sequence] = well
    return whitelist


def parse_well_range(spec):
    """
    Parse a well range specification like 'A1-A3' into a set of wells.
    Currently supports a single range with the same row.
    """
    spec = spec.strip().upper()
    if '-' in spec:
        start, end = spec.split('-')
        row = start[0]
        try:
            start_num = int(start[1:])
            end_num = int(end[1:])
        except ValueError:
            raise ValueError(f"Invalid well numbers in range '{spec}'")
        wells = {f"{row}{i}" for i in range(start_num, end_num + 1)}
    else:
        wells = {spec}
    return wells


def open_output_file(output_dir, sample_id, sample_name, read_label):
    """Open a gzip output file (text mode) for a given sample and read (R1 or R2)."""
    file_path = os.path.join(output_dir, f"${sample_id}_group_{sample_name}_{read_label}.fastq.gz")
    return gzip.open(file_path, 'wt')


def process_fastq(fq1_path, fq2_path, whitelist, barcode_start, barcode_end, output_handles):
    """
    Process paired fastq files.
      - For each read pair, extract a barcode from R2 using the given indices.
      - If the barcode is found in the whitelist, write the read pair to the output.
    """
    stats = {'total': 0, 'matched': 0, 'unmatched': 0}
    with gzip.open(fq1_path, 'rt') as fq1, gzip.open(fq2_path, 'rt') as fq2:
        while True:
            r1_lines = [fq1.readline() for _ in range(4)]
            r2_lines = [fq2.readline() for _ in range(4)]

            # End of file check
            if not r1_lines[1] or not r2_lines[1]:
                break
            stats['total'] += 1

            # Extract barcode from R2 (line 2)
            barcode = r2_lines[1][barcode_start:barcode_end].strip()
            if barcode in whitelist:
                stats['matched'] += 1
                output_handles['R1'].write(''.join(r1_lines))
                output_handles['R2'].write(''.join(r2_lines))
            else:
                stats['unmatched'] += 1
    return stats


def main():
    parser = argparse.ArgumentParser(description="Basic fastq splitter with group filtering")
    parser.add_argument("--sample_id", required=True, help="Sample ID for the output files")
    parser.add_argument("--fq1", required=True, help="Input R1 fastq.gz file")
    parser.add_argument("--fq2", help="Input R2 fastq.gz file; deduced from fq1 if not provided")
    parser.add_argument("--whitelist", required=True, help="CSV file containing whitelist data")
    parser.add_argument("--group", nargs=2, metavar=("SAMPLE", "WELLS"),
                        help="Group specification: sample name and well range (e.g., A1-A3)")
    parser.add_argument("--output", required=True, help="Output directory for split fastq files")
    parser.add_argument("--barcode_start", type=int, default=50, help="Barcode start index (0-based)")
    parser.add_argument("--barcode_end", type=int, default=58, help="Barcode end index (exclusive)")
    args = parser.parse_args()

    # Load the full whitelist from CSV
    whitelist = load_whitelist_csv(args.whitelist)

    # If a group is specified, filter whitelist to only include barcodes from the given wells.
    if args.group:
        sample_name, wells_spec = args.group
        allowed_wells = parse_well_range(wells_spec)
        filtered_whitelist = {barcode: well for barcode, well in whitelist.items() if well in allowed_wells}
        if not filtered_whitelist:
            raise ValueError(f"No barcodes found in whitelist for wells {allowed_wells}")
        whitelist = filtered_whitelist
    else:
        sample_name = "all_samples"

    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)

    # Open output files for the sample (group)
    output_handles = {
        'R1': open_output_file(args.output, args.sample_id, sample_name, "R1"),
        'R2': open_output_file(args.output, args.sample_id, sample_name, "R2")
    }

    # Process the fastq files
    stats = process_fastq(args.fq1, args.fq2, whitelist, args.barcode_start, args.barcode_end, output_handles)

    # Close output files
    output_handles['R1'].close()
    output_handles['R2'].close()

    # Print summary statistics
    print("Processing complete.")
    print(f"Total reads processed: {stats['total']}")
    print(f"Reads matched (within group wells): {stats['matched']}")
    print(f"Reads unmatched: {stats['unmatched']}")

if __name__ == "__main__":
    main()
