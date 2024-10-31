#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=03:00:00
#SBATCH --qos=short
#SBATCH --mem=20G
#SBATCH --job-name starsolo_mapping 


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
seqTech=$3

mapfile -t ACCESSIONS < ${dataDir}/accession_lists/${species}_accessions_v3.txt

# Different STARsolo params depending on sequencing technique
if [[ ${seqTech} == "10xRNAv2" ]]; then
  echo "Mapping using 10xRNAv2 data"
  # cDNA read
  R1="${dataDir}/${species}/fastq/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_3.fastq.gz"
  # Barcode read (cell+UMI)
  R2="${dataDir}/${species}/fastq/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_2.fastq.gz"
  Type="CB_UMI_Simple \
	--soloCBstart 1 \
	--soloCBlen 16 \
	--soloUMIstart 17 \
	--soloUMIlen 10"
  CBwhitelist="${dataDir}/${species}/737K-august-2016.txt" 
  CBmatchWLtype="1MM_multi" # Default
  Strand="Forward" # Default 

elif [[ ${seqTech} == "parse_biosciences" ]]; then
  echo "Mapping using Parse Biosciences data"
  # cDNA read
  R1="${dataDir}/${species}/splitted_fastq_v1/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R1.fastq.gz"
  # Barcode read (cell+UMI)
  R2="${dataDir}/${species}/splitted_fastq_v1/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R2.fastq.gz"
  Type="CB_UMI_Complex \
	--soloCBposition 0_50_0_57 0_30_0_37 0_10_0_17 \
	--soloUMIposition 0_0_0_9"
  CBwhitelist="${dataDir}/${species}/barcodes/bc_data_n26_R1_v3_4 \
	  ${dataDir}/${species}/barcodes/bc_data_v1 \
	  ${dataDir}/${species}/barcodes/bc_data_R3_v3"
  CBmatchWLtype="1MM"
  Strand="Forward" # Default
else
  echo "Sequencing technique unknown: ${seqTech}"
fi

# Creating own subdirectory for results
mkdir -p ${dataDir}/${species}/mapping_splitted_starsolo_v1/results_${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}

# Mapping step and generating count matrix using STAR
STAR --runThreadN 4 \
     --genomeDir ${dataDir}/${species}/genome/genome_index \
     --readFilesCommand zcat \
     --outFileNamePrefix ${dataDir}/${species}/mapping_splitted_starsolo_v1/results_${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}/ \
     --readFilesIn ${R1} ${R2} \
     --soloType ${Type} \
     --soloCBwhitelist ${CBwhitelist} \
     --soloCBmatchWLtype ${CBmatchWLtype} \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand ${Strand} \
     --outSAMattributes CR UR CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within KeepPairs \
     --soloFeatures Gene GeneFull

# Creates the index file
samtools index ${dataDir}/${species}/mapping_splitted_starsolo_v1/results_${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}/*.bam


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds

