#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --qos=vshort
#SBATCH --mem=1G
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
n_cells=$( cat ${mappingDir}/Log.final.out | grep 'input reads' | awk '{print $NF}' )
mapping_rate=$( cat ${mappingDir}/Log.final.out | grep 'Uniquely mapped reads %' | awk '{print $NF}' | sed 's/%//' | awk '{print $1/100}' )
echo ${mapping_rate}

python ${codeDir}/saturation_table.py \
	-b ${bam_file} \
       	-n ${n_cells} \
	-r ${mapping_rate} \
	-o ${output_dir}/output.tsv


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds


