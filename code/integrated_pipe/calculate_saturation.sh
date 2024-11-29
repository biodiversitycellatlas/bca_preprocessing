#!/bin/bash

####################
# define variables #
####################	
mappingDir=$1
codeDir=$2
sample_id=$3

cd ${codeDir}

# prints the folder associated with this run
echo "${mappingDir}"

###############
# run command #
###############
n_cells=$( cat ${mappingDir}/Solo.out/GeneFull/Summary.csv | grep 'Estimated Number of Cells' | sed 's/,/ /g' | awk '{print $NF}' )
n_reads=$( cat ${mappingDir}/Log.final.out | grep 'Number of input reads' | awk '{print $NF}' )
MAPREADS=$( samtools view -F 260 ${mappingDir}/*.bam | wc -l )
map_rate=$( echo "scale=4; ${MAPREADS}/${n_reads}" | bc ) 
temp_folder="${output_dir}/_tmp_${sample_id}"
echo "cells:${n_cells} reads:${n_reads} mapreads:${MAPREADS} maprate:${map_rate}"

python ${codeDir}/saturation_table.py \
	-b ${mappingDir}/*.bam \
       	-n ${n_cells} \
	-r ${map_rate} \
	-t ${temp_folder} \
	-o output.tsv

python ${codeDir}/scripts/plot_curve.py  \
	output.tsv \
	saturation.png \
	--target 0.7 
