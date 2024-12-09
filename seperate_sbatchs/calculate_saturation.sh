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
accession="BCA004_lib_13080AAF_CTTGTAAT-AGTTGGCT"
baseDir="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas"
codeDir="/users/asebe/bvanwaardenburg/git/bca_preprocessing/ext_programs/10x_saturate"
log_out="/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/"

cd ${codeDir}

# prints the accession number associated with this run
echo "${accession}"

mappingDir="${baseDir}/mapping_STARsolo_N/${accession}"
bam_file="${mappingDir}/*.bam"

# create output directory
output_dir="${baseDir}/saturation/${accession}_N"
mkdir -p ${output_dir}


###############
# run command #
###############
n_cells=$( cat ${mappingDir}/*Solo.out/GeneFull/Summary.csv | grep 'Estimated Number of Cells' | sed 's/,/ /g' | awk '{print $NF}' )
n_reads=$( cat ${mappingDir}/*Log.final.out | grep 'Number of input reads' | awk '{print $NF}' )
MAPREADS=$( samtools view -F 260 ${bam_file} | wc -l )
map_rate=$( echo "scale=4; ${MAPREADS}/${n_reads}" | bc ) 
temp_folder="${output_dir}/_tmp_${accession}"
echo "cells:${n_cells} reads:${n_reads} mapreads:${MAPREADS} maprate:${map_rate}"

python ${codeDir}/saturation_table.py \
	-b ${bam_file} \
    -n ${n_cells} \
	-r ${map_rate} \
	-t ${temp_folder} \
	-o ${output_dir}/output.tsv

python ${codeDir}/scripts/plot_curve.py  \
	${output_dir}/output.tsv \
	${output_dir}/saturation.png \
	--target 0.7 

# copy log file to output directory, as it contains the number of reads to reach 0.7 sat
cp ${log_out}*${SLURM_JOB_ID}*.out ${output_dir}

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds


