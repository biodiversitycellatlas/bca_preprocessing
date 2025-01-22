#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=16:00:00
#SBATCH --qos=long
#SBATCH --mem=26G
#SBATCH --job-name cellbender


mapping_files="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF/mapping_STARsolo_N_addattrib/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM/Solo.out/Gene/raw/"
out_dir="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF/CellBender/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM"

# Copy features file as cellbender expects the file to be named genes.tsv
cp ${mapping_files}/features.tsv ${mapping_files}/genes.tsv

# CellBender requires output directory to be created before running
mkdir -p ${out_dir}

cellbender remove-background \
       --input ${mapping_files} \
       --output ${out_dir}/cellbender_output.h5 
