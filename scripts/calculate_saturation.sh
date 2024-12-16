#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --qos=vshort
#SBATCH --mem=10G
#SBATCH --job-name saturation


#################
# # start message #
# #################
# start_epoch=`date +%s`
# echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)


####################
# define variables #
####################	
mappingDir=$1
codeDir=$2
sample_id=$3

# outdir="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002/Saturation/${sample_id}"
# mkdir -p ${outdir}

# prints the folder associated with this run
echo "${mappingDir}"

###############
# run command #
###############
n_cells=$( cat ${mappingDir}/Solo.out/GeneFull/Summary.csv | grep 'Estimated Number of Cells' | sed 's/,/ /g' | awk '{print $NF}' )
n_reads=$( cat ${mappingDir}/Log.final.out | grep 'Number of input reads' | awk '{print $NF}' )
MAPREADS=$( samtools view -F 260 ${mappingDir}/*.bam | wc -l )
map_rate=$( echo "scale=4; ${MAPREADS}/${n_reads}" | bc ) 
temp_folder="_tmp_${sample_id}"
echo "cells:${n_cells} reads:${n_reads} mapreads:${MAPREADS} maprate:${map_rate}"

python ${codeDir}/saturation_table.py \
	-b ${mappingDir}/*.bam \
    -n ${n_cells} \
	-r ${map_rate} \
	-t ${temp_folder} \
	-o output.tsv \
	-d ${codeDir}/scripts

python ${codeDir}/scripts/plot_curve.py  \
	output.tsv \
	saturation.png \
	--target 0.7 