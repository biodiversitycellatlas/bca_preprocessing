#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --qos=vshort
#SBATCH --mem=5G
#SBATCH --job-name saturation


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
accession=$3
codeDir="/users/asebe/bvanwaardenburg/git/10x_saturate"

cd ${codeDir}

# prints the accession number associated with this run
echo "${accession}"

mappingDir="${dataDir}/${species}/mapping_splitted_starsolo/results_${accession}"
bam_file="${mappingDir}/*.bam"

# create output directory
output_dir="${dataDir}/${species}/saturation/${accession}"
mkdir -p ${output_dir}


###############
# run command #
###############
n_cells=$( cat ${mappingDir}/Solo.out/Gene/Summary.csv | grep 'Estimated Number of Cells' | sed 's/,/ /g' | awk '{print $NF}' )
n_reads=$( cat ${mappingDir}/Log.final.out | grep 'Number of input reads' | awk '{print $NF}' )
MAPREADS=$( samtools view -F 260 ${bam_file} | wc -l )
map_rate=$( echo "scale=4; ${MAPREADS}/${n_reads}" | bc ) 

echo "cells:${n_cells} reads:${n_reads} mapreads:${MAPREADS} maprate:${map_rate}"

python ${codeDir}/saturation_table.py \
	-b ${bam_file} \
       	-n ${n_cells} \
	-r ${map_rate} \
	-f ${output_dir}/tags.tab \
	-o ${output_dir}/output.tsv

python ${codeDir}/scripts/plot_curve.py  \
	${output_dir}/output.tsv \
	${output_dir}/saturation.png \
	--target 0.7 

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds


