#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=02:00:00
#SBATCH --qos=short
#SBATCH --mem=20G
#SBATCH --job-name bcl2fastq2


workDir=$1
outDir=$2

# Load bcl2fastq module
module load bcl2fastq2/2.20.0-GCC-13.2.0

# Run bcl2fastq2
/software/sit/EasyBuild/software/bcl2fastq2/2.20.0-GCC-13.2.0/bin/bcl2fastq \
  --runfolder-dir ${workDir} \
  --output-dir ${outDir} \
  --no-lane-splitting \
  --create-fastq-for-index-reads
