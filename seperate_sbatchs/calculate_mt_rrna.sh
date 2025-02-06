#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --qos=vshort
#SBATCH --mem=6G
#SBATCH --job-name featureCounts

# ------------------------------------------------------------------
# Define Variables
# ------------------------------------------------------------------
# e.g. work_dir="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas_sep/"  
work_dir=$1
annot_type="GTF"
ref_star_gtf="${work_dir}/genome/genomic.gtf"

output_dir="${work_dir}/FeatureCounts/"
mkdir -p ${output_dir}


# ------------------------------------------------------------------
# Check for GFF/GTF
# Editing the GFF file format if needed
# ------------------------------------------------------------------
if [ "${annot_type}" == "GFF" ]; then
    echo "Checking whether the GFF file is correctly formatted..."
    cp ${ref_star_gtf} ${output_dir}/
    cat ${output_dir}/*.gff | sed 's/ID=/gene_id=/g' > ${output_dir}/genomic_v2.gff
    new_ref="${output_dir}/genomic_v2.gff"
else
    new_ref=${ref_star_gtf}
fi


for config_dir in "${work_dir}"/mapping_STARsolo/*; do
    [[ ! -d "${config_dir}" ]] && continue  # Skip if not a directory
    for data_dir in "${config_dir}"/*; do
        [[ ! -d "${data_dir}" ]] && continue  # Skip if not a directory

        config="${config_dir##*_}"
        data_name=$( basename ${data_dir} )
        echo "data_dir: ${data_dir}, ${data_name}, config_dir: ${config_dir}, ${config}"
        
        output_dir="${work_dir}/FeatureCounts/FeatureCounts_${config}/${data_name}"
        mkdir -p ${output_dir}

        # Create seperate BAM with filtered results - filter for multimapped reads only present in primary alignment
        samtools view -h -F 256 ${data_dir}/Aligned.sortedByCoord.out.bam | grep -E "^\@|NH:i:[2-9]" | samtools view -b -o ${output_dir}/multimapped_primealign.bam
        # Create seperate BAM with filtered results - filter for multimapped reads present in both secondary and primary alignments
        samtools view -h ${data_dir}/Aligned.sortedByCoord.out.bam | grep -E "^\@|NH:i:[2-9]" | samtools view -b -o ${output_dir}/multimapped_allalign.bam

        # Calculation rRNA
        featureCounts -t rRNA -a ${new_ref} -o ${output_dir}/feat_counts_rRNA.txt ${data_dir}/Aligned.sortedByCoord.out.bam
        featureCounts -M --fraction -t rRNA -a ${new_ref} -o ${output_dir}/feat_counts_rRNA_mmpa.txt ${output_dir}/multimapped_primealign.bam
        featureCounts -M --fraction -t rRNA -a ${new_ref} -o ${output_dir}/feat_counts_rRNA_mmaa.txt ${output_dir}/multimapped_allalign.bam

        # Calculation mtDNA
        grep 'mtDNA' ${new_ref} > ${output_dir}/mtDNA_only.gtf
        featureCounts -t exon -a ${output_dir}/mtDNA_only.gtf -o ${output_dir}/feat_counts_mtDNA.txt ${data_dir}/Aligned.sortedByCoord.out.bam
        featureCounts -M --fraction -t exon -a ${output_dir}/mtDNA_only.gtf -o ${output_dir}/feat_counts_mtDNA_mmpa.txt ${output_dir}/multimapped_primealign.bam
        featureCounts -M --fraction -t exon -a ${output_dir}/mtDNA_only.gtf -o ${output_dir}/feat_counts_mtDNA_mmaa.txt ${output_dir}/multimapped_allalign.bam

        # Per-gene counts for multimappers
        featureCounts -M --fraction -a ${new_ref} -o ${output_dir}/feat_counts_multimappers_genelevel.txt ${output_dir}/multimapped_allalign.bam

        # Saves feature counts as clean gene counts table
        tail -n +3 ${output_dir}/feat_counts_multimappers_genelevel.txt | cut -f1,7 > ${output_dir}/multimapper_gene_counts.txt

        # Sort the counts
        sort -k2,2nr ${output_dir}/multimapper_gene_counts.txt > ${output_dir}/sorted_multimapper_gene_counts.txt

    done
done