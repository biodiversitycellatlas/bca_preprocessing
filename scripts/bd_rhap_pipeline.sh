#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=22:00:00
#SBATCH --qos=long
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name bdrhap_pipeline_BCA005_NvecV6
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=vanwaardenburgbonita@gmail.com

####################
# define variables #
####################
# res_dir=$1
# sample_id=$2

# echo "read 1: ${read_1}"
# echo "read 2: ${read_2}"
# echo "barcode file: ${barcode_path}"


workDir="/users/asebe/bvanwaardenburg/CRSwDev-cwl-8d252728f875/v2.2.1"
outDir="/users/asebe/bvanwaardenburg/git/data/241106_BD_Rhapsody_Nvec/Nvec_BDpipeline_seperate/BCA005_NvecV6"

cd ${workDir}

cwl-runner \
    --outdir ${outDir} \
    --singularity \
    rhapsody_pipeline_2.2.1.cwl \
    pipeline_inputs_BCA005_2.2.1.yml