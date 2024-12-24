#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=00:40:00
#SBATCH --qos=vshort
#SBATCH --mem=50G
#SBATCH --job-name souporcell 

resDir="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF/"
bc_file="${resDir}/mapping_STARsolo_N/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM/Solo.out/GeneFull/raw/barcodes.tsv"
bam_file="${resDir}/mapping_STARsolo_N/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM/Aligned.sortedByCoord.out.bam"
ref_fasta="${resDir}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta"

baseDir="/users/asebe/bvanwaardenburg/git/bca_preprocessing"

echo "Barcodes file: \${bc_file}"
echo "BAM file: \${bam_file}"

${baseDir}/ext_programs/souporcell/souporcell_pipeline.py \
    --bam ${bam_file} \
    --barcodes ${bc_file} \
    --fasta ${ref_fasta} \
    --threads 40 \
    --out_dir ${resDir}/doublet_det \
    --clusters 4