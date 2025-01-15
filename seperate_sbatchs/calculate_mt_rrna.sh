#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=00:20:00
#SBATCH --qos=vshort
#SBATCH --mem=5G
#SBATCH --job-name featureCounts 


#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

###############
# run command #
###############

resDir="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF"
gtf_file="${resDir}/genome/Nvec_v5_merged_annotation_sort.gtf"


# For each STARsolo mapping directory
for map_dir in ${resDir}/mapping_STARsolo_N_addattrib; do
    echo "mapping dir: ${map_dir}"

    # For each sample directory under the mapping directory
    for sample_dir in $map_dir/*; do
      echo "sample dir: ${sample_dir} "
      
      # Create seperate BAM with filtered results - filter for multimapped reads only present in primary alignment
      samtools view -h -F 256 ${sample_dir}/Aligned.sortedByCoord.out.bam | grep -E "^\@|NH:i:[2-9]" | samtools view -b -o ${sample_dir}/multimapped_primealign.bam
      # Create seperate BAM with filtered results - filter for multimapped reads present in both secondary and primary alignments
      samtools view -h ${sample_dir}/Aligned.sortedByCoord.out.bam | grep -E "^\@|NH:i:[2-9]" | samtools view -b -o ${sample_dir}/multimapped_allalign.bam

      # rRNA
      featureCounts -t rRNA -a ${gtf_file} -o ${sample_dir}/feat_counts_rRNA.txt ${sample_dir}/Aligned.sortedByCoord.out.bam
      featureCounts -M --fraction -t rRNA -a ${gtf_file} -o ${sample_dir}/feat_counts_rRNA_mmpa.txt ${sample_dir}/multimapped_primealign.bam
      featureCounts -M --fraction -t rRNA -a ${gtf_file} -o ${sample_dir}/feat_counts_rRNA_mmaa.txt ${sample_dir}/multimapped_allalign.bam
      
      # mtDNA
      grep 'mtDNA' ${gtf_file} > ${sample_dir}/mtDNA_only.gtf
      featureCounts -t exon -a ${sample_dir}/mtDNA_only.gtf -o ${sample_dir}/feat_counts_mtDNA.txt ${sample_dir}/Aligned.sortedByCoord.out.bam
      featureCounts -M --fraction -t exon -a ${sample_dir}/mtDNA_only.gtf -o ${sample_dir}/feat_counts_mtDNA_mmpa.txt ${sample_dir}/multimapped_primealign.bam
      featureCounts -M --fraction -t exon -a ${sample_dir}/mtDNA_only.gtf -o ${sample_dir}/feat_counts_mtDNA_mmaa.txt ${sample_dir}/multimapped_allalign.bam

      # Per-gene counts for multimappers
      featureCounts -M --fraction -a ${gtf_file} -o ${sample_dir}/feat_counts_multimappers_genelevel.txt ${sample_dir}/multimapped_allalign.bam
      # Saves feature counts as clean gene counts table
      tail -n +3 ${sample_dir}/feat_counts_multimappers_genelevel.txt | cut -f1,7 > ${sample_dir}/multimapper_gene_counts.txt
      # Sort the counts
      sort -k2,2nr ${sample_dir}/multimapper_gene_counts.txt > ${sample_dir}/sorted_multimapper_gene_counts.txt

    done
done


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds

