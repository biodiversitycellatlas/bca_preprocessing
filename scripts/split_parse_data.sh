#!/bin/bash

####################
# define variables #
####################
res_dir=$1
sample_id=$2
read_1=$3
read_2=$4
barcode_path=$5		

echo "read 1: ${read_1}"
echo "read 2: ${read_2}"
echo "barcode file: ${barcode_path}"

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
  
  barcodes_r1=$( cat ${barcode_path} | grep -E ",A[${num_wells}]," | awk -F ',' '{print $2}' | sed ':a;N;$!ba;s/\n/\|/g' )
  echo "barcodes matched: ${barcodes_r1}"

  # Step 1: Select sequences matching the barcodes and write to file
  # version 2 allows for n number of mismatces in the barcode
  echo "${read_2}"
  zcat "${read_2}" | grep --no-group-separator -B 1 -A 2 -E "${barcodes_r1}" > "${sample_id}_${name}_R2.fastq"
  echo "Step 1 completed, created ${sample_id}_${name}_R2.fastq"

  # Step 2: Extract headers from the filtered R2 FASTQ file
  cat "${sample_id}_${name}_R2.fastq" | grep '^@' |  sed 's/^@//g' | awk '{print $1}' > "${sample_id}_${name}_headers_R2.txt"
  echo "Step 2 completed, created filtered_headers_R2.txt"

  # Step 3
  seqtk subseq ${read_1} "${sample_id}_${name}_headers_R2.txt" > "${sample_id}_${name}_R1.fastq"
  echo "Step 3 completed, created ${sample_id}_${name}_R1.fastq"

  # Step 4: Remove intermediate files
  rm "${sample_id}_${name}_headers_R2.txt"
  gzip "${sample_id}_${name}_R2.fastq"
  gzip "${sample_id}_${name}_R1.fastq"

done 6< ${sw_file}


