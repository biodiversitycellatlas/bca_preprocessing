#!/bin/bash

####################
# define variables #
####################
res_dir=$1
sample_id=$2
fastq_files=$3
barcode_path=$4		# directory with the barcode files for each round
output_dir=$5

echo "============================"
echo "res_dir: ${res_dir}"
echo "sample_id: ${sample_id}"
echo "fastq_files: ${fastq_files}"
echo "barcode_path: ${barcode_path}"
echo "output_dir: ${output_dir}"
echo "============================"

# retrieve FASTQ file names from array
read_1=${fastq_files[0]}
read_2=${fastq_files[1]}

# selects the barcode files (for chemistry v3!!!)
bc_1="${barcode_path}/bc_data_n26_R1_v3_4.csv"
bc_2="${barcode_path}/bc_data_v1.csv"
bc_3="${barcode_path}/bc_data_R3_v3.csv"

# checks which sample_wells file to use depending on the name 
sample_num=$( echo ${sample_id} | grep -oP '(?<=BCA00)\d' | head -1 )
echo "sample num: ${sample_num}, resdir: ${res_dir}"
if [ "${sample_num}" -lt 3 ];
then
  sw_file="${res_dir}/B001_B002_samplewells.txt"
else
  sw_file="${res_dir}/B003_B004_samplewells.txt"
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
  echo "${name}: ${num_wells}"
  
  barcodes_r1=$( cat ${bc_1} | grep -E ",A[${num_wells}]," | awk -F ',' '{print $2}' | sed ':a;N;$!ba;s/\n/\|/g' )
  # barcodes_r2=$( cat $bc_2 | grep -E ",A[${num_wells}]," | awk -F',' '{print $2}' | sed ':a;N;$!ba;s/\n/\|/g' )
  # barcodes_r3=$( cat $bc_3 | grep -E ",A[${num_wells}]," | awk -F',' '{print $2}' | sed ':a;N;$!ba;s/\n/\|/g' )
  # echo "round 1: '${barcodes_r1}' \n round 2: '${barcodes_r2}' \n round 3: '${barcodes_r3}'"
  
  # Step 1: Select sequences matching the barcodes and write to file
  # version 2 allows for n number of mismatces in the barcode
  echo "${read_2}"
  zcat "${read_2}" | grep --no-group-separator -B 1 -A 2 -E "${barcodes_r1}"| gzip  > "${output_dir}/${sample_id}_${name}_R2.fastq.gz"
  echo "Step 1 completed, created ${output_dir}/${sample_id}_${name}_R2.fastq.gz"

  # Step 2: Extract headers from the filtered R2 FASTQ file
  zcat "${output_dir}/${sample_id}_${name}_R2.fastq.gz" | grep '^@' |  sed 's/^@//g' | awk '{print $1}' > "${output_dir}/${sample_id}_${name}_headers_R2.txt"
  echo "Step 2 completed, created ${output_dir}/filtered_headers_R2.txt"

  # Step 3
  seqtk subseq ${read_1} "${output_dir}/${sample_id}_${name}_headers_R2.txt" | gzip > "${output_dir}/${sample_id}_${name}_R1.fastq.gz"
  echo "Step 3 completed, created ${output_dir}/${sample_id}_${name}_R1.fastq.gz"

  # Step 4: Remove intermediate files
  rm "${output_dir}/${sample_id}_${name}_headers_R2.txt"

done 6< ${sw_file}


