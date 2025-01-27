#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=10:00:00
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --job-name kraken2_db

# Define variables
dataDir="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas_sep/mapping_STARsolo_N/BCA004_lib_13080AAF_CTTGTAAT-AGTTGGCT_all"

# Build Kraken2 default database 11h + 26G
kraken2-build --standard \
    --db /users/asebe/bvanwaardenburg/git/data/kraken2_dir/standard_db/ \
    --threads 32 \
    --max-db-size 8000000000

# Randomly subsample 10% of the data 30 min + 26G
# seqkit sample -p 0.1 \
#     -o ${dataDir}/subsampled_unmapped_reads.Trinity.fasta \
#     ${dataDir}/Trinity_Assembled_unmapped_reads.Trinity.fasta

# Perform taxonomic classification
# kraken2 --db /users/asebe/bvanwaardenburg/git/data/kraken2_dir/standard_db/ \
#    --threads 8 \
#    --output /users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas_sep/kraken2/results.txt \
#    --report /users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas_sep/kraken2/report.txt \
#    ${dataDir}/subsampled_unmapped_reads.Trinity.fasta
