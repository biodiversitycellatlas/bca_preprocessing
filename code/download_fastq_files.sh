#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --array=1-28
#SBATCH --time=01:30:00
#SBATCH --qos=shorter
#SBATCH --mem=8G
#SBATCH --job-name download_fastq 


#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)


##################################
# make bash behave more robustly #
##################################
set -e
set -u
set -o pipefail


###################
# set environment #
###################
module load SRA-Toolkit/3.1.1-gompi-2023b


###############
# run command #
###############
mapfile -t ACCESSIONS < accession_list.txt
file_path="/users/asebe/bvanwaardenburg/git/bca_preprocessing/data/spol/fastq"
fastq-dump --split-files --gzip --outdir "$file_path" ${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

