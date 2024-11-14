#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=02:00:00
#SBATCH --qos=shorter
#SBATCH --mem=16G
#SBATCH --job-name split_demult 


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
codeDir="/users/asebe/bvanwaardenburg/git/SPLiT-Seq_demultiplexing"
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
output_dir="${dataDir}/${species}/splitted_fastq_demult"
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

bash ${codeDir}/splitseqdemultiplex_0.2.2.sh \
	-v split \
	--round1barcodes ${bc_1} \
	--round2barcodes ${bc_2} \
	--round3barcodes ${bc_3} \
	--fastqF ${read_1} \
	--fastqR ${read_2} \
        --outputdir ${output_dir} \
	--targetMemory 16000 


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds


