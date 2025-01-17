#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=3:00:00
#SBATCH --qos=short
#SBATCH --mem=26G
#SBATCH --job-name souporcell


mapping_files="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF/mapping_STARsolo_N_addattrib/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM"
bc_file=$(ls ${mapping_files}/Solo.out/Gene/raw/barcodes.tsv | head -n 1)
bam_file=$(ls ${mapping_files}/Aligned.sortedByCoord.out.bam | head -n 1)
ref_fasta="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF/genome/Nvec_vc1.1_gDNA_mtDNA.fasta" 
out_dir="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF/souporcell"

echo "Barcodes file: ${bc_file}"
echo "BAM file: ${bam_file}"
    
singularity exec /users/asebe/bvanwaardenburg/souporcell_release.sif souporcell_pipeline.py \
	--bam ${bam_file} \
        --barcodes ${bc_file} \
        --fasta ${ref_fasta} \
        --threads 40 \
        --out_dir ${out_dir} \
	--clusters 12
