#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=02:30:00
#SBATCH --qos=short
#SBATCH --mem=26G
#SBATCH --job-name geneext

sample="AACGAACTGT"
out_dir="/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/geneExt/${sample}"

# Create output directory if it does not exist
mkdir -p ${out_dir}
cd ${out_dir}

# Run geneext
python /users/asebe/bvanwaardenburg/git/bca_preprocessing/ext_programs/GeneExt/geneext.py \
    -g ../../genome/Nvec_v6_merged_noorphanpeaks_sort.gtf \
    -b ../../mapping_STARsolo/mapping_STARsolo_fastqsplitted_GTFv6/${sample}/${sample}_Aligned.sortedByCoord.out.bam \
    -o ${sample}_geneext.gtf \
    -j 4
