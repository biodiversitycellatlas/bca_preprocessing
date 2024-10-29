#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --qos=vshort
#SBATCH --mem=20G
#SBATCH --job-name starsolo_index 


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

# Retrieve the first accession number
first_accs=$(head -1 ${dataDir}/accession_lists/${species}_accessions.txt)

# Different file-ending depending on sequencing technique
if [[ "${seqTech}" == "10xRNAv2" ]]; then
  first_fastq="${dataDir}/${species}/fastq/${first_accs}_3.fastq.gz" 
else
  first_fastq="${dataDir}/${species}/fastq/${first_accs}_R1_001.fastq.gz"  
fi

# Calculate SJBD overhang using the first read from the first fastq file
sjdb_overhang=$(zcat "${first_fastq}" 2>/dev/null | awk 'NR==2 {print length($0)-1; exit}' || echo "") 

# Create directory where genome index will be stored
mkdir -p ${dataDir}/${species}/genome/genome_index

# Specify file paths
genomeFastaFile=$(ls ${dataDir}/${species}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta)
GTFfile=$(ls ${dataDir}/${species}/genome/Nvec_v5_merged_annotation_sort.gtf) 

# Generating genome index using STAR
STAR \
  --runMode genomeGenerate \
  --genomeDir "${dataDir}/${species}/genome/genome_index" \
  --genomeFastaFiles "${genomeFastaFile}" \
  --sjdbGTFfile "${GTFfile}" \
  --sjdbOverhang "${sjdb_overhang}" \
  --genomeSAindexNbases 12

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds

