#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=11:00:00
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --job-name demux_umitools

cd /users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec

umi_tools extract \
    --extract-method='string' \
    --bc-pattern='NNNNNNNNNN' \
    --stdin=bcl2fastq2_nolanesplit_index/Undetermined_S0_I2_001.fastq.gz \
    --read2-in=bcl2fastq2_nolanesplit_index/Undetermined_S0_R2_001.fastq.gz  \
    --read2-out=demux_UMItools/R2_i5_extracted.fastq.gz \
    --log demux_UMItools/umi_extract_i5.log


# umi_tools extract \
#     --extract-method='regex' \
#     --bc-pattern='(?P<cell_1>.{16})(?P<umi_1>.{10})' \
#     --stdin=bcl2fastq2_nolanesplit_index/Undetermined_S0_R1_001.fastq.gz \
#     --read2-in=bcl2fastq2_nolanesplit_index/Undetermined_S0_R2_001.fastq.gz  \
#     --read2-out=demux_UMItools/R2_demux_final.fastq.gz \
#     --log demux_UMItools/umi_extract_R1.log

