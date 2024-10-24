#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=01:00:00
#SBATCH --qos=shorter
#SBATCH --mem=1G
#SBATCH --job-name split_fastqs 


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
barcode_path=$3		# directory with the barcode files for each round

# reads the accession list into memory
mapfile -t ACCESSIONS < ${dataDir}/accession_lists/${species}_accessions.txt

# prints the accession number associated with this run
echo "${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}"

# selects the fastq files
R1="${dataDir}/${species}/fastq/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R1_001.fastq.gz"
R2="${dataDir}/${species}/fastq/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R2_001.fastq.gz"

# selects the barcode files (for chemistry v3!!!)
bc_1="${barcode_path}/bc_data_n26_R1_v3_4.csv"
bc_2="${barcode_path}/bc_data_v1.csv"
bc_3="${barcode_path}/bc_data_R3_v3.csv"

# create output directory
output_dir="${dataDir}/${species}/splitted_fastq"
mkdir -p ${output_dir}

# checks which sample_wells file to use depending on the run 
if [ "$SLURM_ARRAY_TASK_ID" -lt 3 ];
then
  sw_file="${dataDir}/${species}/B001_B002_samplewells.txt"
else
  sw_file="${dataDir}/${species}/B003_B004_samplewells.txt"
fi

###############
# run command #
###############
# for x in sw_file (e.g.: 'nvec_DSP A1-A3')
	# count how many matches were found for each round, report it in a file and say how many reads were not assigned
	# name=x[0]
	# wells_alpha=x[1].alpha()   	split letters and numbers
	# wells_num=x[1].num()  	e.g. 1-3
	# barcodes_r1=$(cat $bc_1 | grep -E ",A[$wells_num]," | awk -F',' '{print $2}')
	# barcodes_r2=$(cat $bc_2 | grep -E ",A[$wells_num]," | awk -F',' '{print $2}')
	# barcodes_r3=$(cat $bc_3 | grep -E ",A[$wells_num]," | awk -F',' '{print $2}')
	# files=$(zcat $R1, $R2 | grep -E $barcodes_r1 | grep -E $barcodes_r2 | grep -E $barcodes_r3 > ${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R1_${name}.fastq.gz, ${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R1_${name}.fastq.gz )	


# read counters
total_reads=0
unassigned_reads=0

# loops through the sample_wells file (e.g. line: 'DSP A1-A3'), -u 6 set as unknown file descriptor 
## TODO: check if the letters are the same, and change it from A to dynamic
while read -u 6 line;
do
  name=$( echo "${line}" | awk '{print $1}' )
  wells=$( echo "${line}" | awk '{print $2}' )
  alpha_wells=$( echo ${wells} | tr '-' ' ' | sed 's/[0-9]//g' )
  num_wells=$( echo ${wells} | sed 's/[A-Za-z]//g' )
  echo "${line}: ${name} ${wells} num_wells : ${num_wells} alpha_wells:${alpha_wells}"
  
  # change sed argument
  barcodes_r1=$( cat $bc_1 | grep -E ",A[${num_wells}]," | awk -F',' '{print $2}' | sed ':a;N;$!ba;s/\n/\/|/g' )
  barcodes_r2=$( cat $bc_2 | grep -E ",A[${num_wells}]," | awk -F',' '{print $2}' | sed ':a;N;$!ba;s/\n/\/|/g' )
  barcodes_r3=$( cat $bc_3 | grep -E ",A[${num_wells}]," | awk -F',' '{print $2}' | sed ':a;N;$!ba;s/\n/\/|/g' )
  echo "round 1: ${barcodes_r1} \n round 2: ${barcodes_r2} \n round 3: ${barcodes_r3}"
  
  # select the sequences that match the combination of the 3 barcodes, and write to file
  filtered_R2=$( zcat ${R2} | grep -E '${barcodes_r1}' | grep -E '${barcodes_r2}' | grep -E '${barcodes_r3}' > ${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R2_${name}.fastq )
 
done 6< ${sw_file}

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds


