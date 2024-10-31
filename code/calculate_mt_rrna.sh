#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=00:20:00
#SBATCH --qos=vshort
#SBATCH --mem=5G
#SBATCH --job-name featureCounts 


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
accession=$3

output_dir="${dataDir}/${species}/mapping_splitted_starsolo_v1/results_${accession}"
gtf_file="${dataDir}/${species}/genome/Nvec_v5_merged_annotation_sort.gtf"

featureCounts -M --fraction -f \
	-t rRNA \
	-a ${gtf_file} \
	-o ${output_dir}/feat_counts_rRNA.txt \
	${output_dir}/Aligned.sortedByCoord.out.bam


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds

