#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=../logs/%x.%j.out
#SBATCH --error=../logs/%x.%j.err
#SBATCH --time=01:00:00
#SBATCH --qos=vshort
#SBATCH --mem=32G
#SBATCH --job-name bcl_oakseq

output_dir=$1

# Load bcl2fastq module
module load bcl2fastq2/2.20.0-GCC-13.2.0

# Run bcl2fastq2
# --no-lane-splitting:                Do not split the output by lane
# --create-fastq-for-index-reads:     Create FASTQ files for index reads
/software/sit/EasyBuild/software/bcl2fastq2/2.20.0-GCC-13.2.0/bin/bcl2fastq \
  --runfolder-dir ${output_dir}/raw_data/ \
  --output-dir ${output_dir}/fastq/ \
  --no-lane-splitting \
  --create-fastq-for-index-reads