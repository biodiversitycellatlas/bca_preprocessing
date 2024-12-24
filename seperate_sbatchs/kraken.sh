#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=02:30:00
#SBATCH --qos=shorter
#SBATCH --mem=50G
#SBATCH --job-name kraken2 

kraken2-build --standard --threads 32 --db /users/asebe/bvanwaardenburg/git/bca_preprocessing/ext_programs/kraken2_dir/standard_db/

kraken2 \
    -db /users/asebe/bvanwaardenburg/git/bca_preprocessing/ext_programs/kraken2_dir/standard_db/ \
    /users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas_sep/mapping_STARsolo_N/BCA004_lib_13080AAF_CTTGTAAT-AGTTGGCT_all/Trinity_Assembled_unmapped_reads.Trinity.fasta