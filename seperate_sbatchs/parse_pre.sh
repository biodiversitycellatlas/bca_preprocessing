#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=10:00:00
#SBATCH --qos=normal
#SBATCH --mem=30G
#SBATCH --job-name parse_pre 


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


####################
# define variables #
####################
species=$1
dataDir=$2

mapfile -t ACCESSIONS < ${dataDir}/accession_lists/${species}_accessions.txt
echo "${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}"

mode="pre"
R1="${dataDir}/${species}/fastq/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R1_001.fastq.gz"
R2="${dataDir}/${species}/fastq/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R2_001.fastq.gz"
chem="v3"

if [ $SLURM_ARRAY_TASK_ID -lt 3 ];
then
  sample_file="${dataDir}/${species}/B001_B002_samplewells.txt"
else
  sample_file="${dataDir}/${species}/B003_B004_samplewells.txt"
fi

outdir="${dataDir}/${species}/corrected_fastq_long/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}"
ref_genome="${dataDir}/${species}/genome/parse_refgenome"


###############
# run command #
###############
echo "running Parse Biosciences pipeline for library ${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}"
split-pipe -m ${mode} --chemistry ${chem} --fq1 ${R1} --fq2 ${R2} \
	--output_dir ${outdir} --genome_dir ${ref_genome} \
	--samp_list ${sample_file} 


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds

