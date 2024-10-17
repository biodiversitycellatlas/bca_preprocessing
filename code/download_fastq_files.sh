#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=02:30:00
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


###############
# run command #
###############
species=$1
dataDir=$2

mapfile -t ACCESSIONS < ${dataDir}/accession_lists/${species}_accessions.txt
fastq-dump --split-files --gzip --outdir ${dataDir}/${species}/fastq \
	${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds

