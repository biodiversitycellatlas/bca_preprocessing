#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=01:30:00
#SBATCH --qos=shorter
#SBATCH --mem=30G
#SBATCH --job-name split_fastqs1 


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
read_1="${dataDir}/${species}/fastq/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R1_001.fastq.gz"
read_2="${dataDir}/${species}/fastq/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_R2_001.fastq.gz"

# selects the barcode files (for chemistry v3!!!)
bc_1="${barcode_path}/bc_data_n26_R1_v3_4.csv"
bc_2="${barcode_path}/bc_data_v1.csv"
bc_3="${barcode_path}/bc_data_R3_v3.csv"

# create output directory
output_dir="${dataDir}/${species}/splitted_fastq_v1"
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

# loops through the sample_wells file (e.g. line: 'DSP A1-A3'), -u 6 set as unknown file descriptor 
## TODO: check if the letters are the same, and change it from A to dynamic
while read -u 6 line;
do
  name=$( echo "${line}" | awk '{print $1}' )
  wells=$( echo "${line}" | awk '{print $2}' )
  alpha_wells=$( echo ${wells} | tr '-' ' ' | sed 's/[0-9]//g' )
  num_wells=$( echo ${wells} | sed 's/[A-Za-z]//g' )
  echo "${line}: ${name} ${wells} num_wells : ${num_wells} alpha_wells:${alpha_wells}"
  
  barcodes_r1=$( cat $bc_1 | grep -E ",A[${num_wells}]," | awk -F',' '{print $2}' | sed ':a;N;$!ba;s/\n/\|/g' )
  barcodes_r2=$( cat $bc_2 | grep -E ",A[${num_wells}]," | awk -F',' '{print $2}' | sed ':a;N;$!ba;s/\n/\|/g' )
  barcodes_r3=$( cat $bc_3 | grep -E ",A[${num_wells}]," | awk -F',' '{print $2}' | sed ':a;N;$!ba;s/\n/\|/g' )
  echo "round 1: '${barcodes_r1}' \n round 2: '${barcodes_r2}' \n round 3: '${barcodes_r3}'"
  
  # Step 1: Select sequences matching the barcodes and write to file
  # version 2 allows for n number of mismatces in the barcode
  echo "${read_2}"
  zcat "${read_2}" | grep -B 1 -A 2 -E "${barcodes_r1}"| gzip  > "${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_${name}_R2.fastq.gz"
  echo "Step 1 completed, created ${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_${name}_R2.fastq.gz"

  # Step 2: Extract headers from the filtered R2 FASTQ file
  zcat "${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_${name}_R2.fastq.gz" | grep '^@' |  sed 's/^@//g' | awk '{print $1}' > "${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_${name}_headers_R2.txt"
  echo "Step 2 completed, created ${output_dir}/filtered_headers_R2.txt"

  # Step 3
  seqtk subseq ${read_1} "${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_${name}_headers_R2.txt" | gzip > "${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_${name}_R1.fastq.gz"
  echo "Step 3 completed, created ${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_${name}_R1.fastq.gz"

  # Step 4: Remove intermediate files
  rm "${output_dir}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}_${name}_headers_R2.txt"

done 6< ${sw_file}

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds


