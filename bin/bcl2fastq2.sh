#!/bin/bash

workDir=$1
outDir=$2

# Load bcl2fastq module
module load bcl2fastq2/2.20.0-GCC-13.2.0

# Run bcl2fastq2
# --no-lane-splitting: Do not split the output by lane
# --create-fastq-for-index-reads: Create FASTQ files for index reads
/software/sit/EasyBuild/software/bcl2fastq2/2.20.0-GCC-13.2.0/bin/bcl2fastq \
  --runfolder-dir ${workDir} \
  --output-dir ${outDir} \
  --no-lane-splitting \
  --create-fastq-for-index-reads
