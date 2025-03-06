#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=2:00:00
#SBATCH --qos=shorter
#SBATCH --mem=32G
#SBATCH --job-name demux_oakseq

####################
# define variables #
####################
# res_dir=$1
# sample_id=$2

# echo "read 1: ${read_1}"
# echo "read 2: ${read_2}"
# echo "barcode file: ${barcode_path}"


read_1="/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/bcl2fastq2_nolanesplit_index/Undetermined_S0_R1_001.fastq.gz"
read_i5="/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/bcl2fastq2_nolanesplit_index/Undetermined_S0_I2_001.fastq.gz"
out_file="/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/demux_UMItools/combined_R1_I2.fastq.gz"

python3 demux_oakseq.py \
    ${read_1} \
    ${read_i5} \
    ${out_file}