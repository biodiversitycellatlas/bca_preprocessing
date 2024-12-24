#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=09:00:00
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --job-name trinity

module load GCC

Trinity \
  --seqType fa \
  --max_memory 500G \
  --single unmapped.fasta \
  --CPU 32 \
  --SS_lib_type F \
  --output Trinity_Assembled_unmapped_reads \
  --full_cleanup

