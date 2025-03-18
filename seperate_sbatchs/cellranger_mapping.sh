#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=05:00:00
#SBATCH --qos=short
#SBATCH --mem=26G
#SBATCH --job-name cellranger


# Add CellRanger to PATH
export PATH=/users/asebe/bvanwaardenburg/cellranger-9.0.1:$PATH

out_dir="/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/mapping_CellRanger/Nvec_reference_v4_extended"

cd ${out_dir}

cellranger count \
    --id=TTCGACAAGC_count \
    --transcriptome=../../genome/cellranger_ref/Nvec_reference_v4_extended/ \
    --fastqs=../../fastq_splitted/ \
    --sample=TTCGACAAGC \
    --chemistry=auto \
    --create-bam true