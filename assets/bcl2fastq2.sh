#!/bin/bash

# Assumes bcl2fastq2 module is available, e.g.:
# module load bcl2fastq2

# Input arguments
input_dir=$1
output_dir=$2

# Create a fake sample sheet to add p5 + p7 indexes to the headers
echo -e "[DATA]\nLane,Sample_ID,Sample_Name,index,index2\n,fake,fake,NNNNNNNNNN,NNNNNNNNNN" > fake_sample_sheet.csv

# Run bcl2fastq2 with the fake sample sheet
bcl2fastq \
  --runfolder-dir ${input_dir} \
  --output-dir ${output_dir} \
  --sample-sheet fake_sample_sheet.csv \
  --barcode-mismatches 1 --ignore-missing-positions --ignore-missing-controls --ignore-missing-filter \
  --ignore-missing-bcls --no-lane-splitting --minimum-trimmed-read-length 15 --mask-short-adapter-reads 15
