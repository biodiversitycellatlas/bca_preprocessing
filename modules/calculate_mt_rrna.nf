// =============  CALCULATION rRNA & mtRNA  ============= \\ 
// 

process CALC_MT_RRNA {   
    publishDir "${params.resDir}/rRNA_mtDNA/${sample_id}_${config_name}", mode: 'copy'
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    output:
    path("*")     
       
    script:
    """
    echo "\n\n==================  CALCULATION rRNA & mtDNA ${config_name} =================="
    echo "Sample ID: ${sample_id}"
    echo "Config name: ${config_name}"
    echo "Mapping files: ${mapping_files}"

    # Create seperate BAM with filtered results - filter for multimapped reads only present in primary alignment
    samtools view -F 256 Aligned.sortedByCoord.out.bam | grep -E "^\@|NH:i:[2-9]" > multimapped_primealign.bam
    # Create seperate BAM with filtered results - filter for multimapped reads present in both secondary and primary alignments
    samtools view Aligned.sortedByCoord.out.bam | grep -E "^\@|NH:i:[2-9]" > multimapped_allalign.bam

    # Calculation rRNA
    featureCounts -t rRNA -a ${params.ref_star_gtf} -o feat_counts_rRNA.txt Aligned.sortedByCoord.out.bam
    featureCounts -t rRNA -a ${params.ref_star_gtf} -o feat_counts_rRNA_mmpa.txt multimapped_primealign.bam
    featureCounts -t rRNA -a ${params.ref_star_gtf} -o feat_counts_rRNA_mmaa.txt multimapped_allalign.bam
    
    # Calculation mtDNA
    grep 'mtDNA' ${params.ref_star_gtf} > mtDNA_only.gtf
    featureCounts -t exon -a mtDNA_only.gtf -o feat_counts_mtDNA.txt Aligned.sortedByCoord.out.bam
    featureCounts -t exon -a mtDNA_only.gtf -o feat_counts_mtDNA_mmpa.txt multimapped_primealign.bam
    featureCounts -t exon -a mtDNA_only.gtf -o feat_counts_mtDNA_mmaa.txt multimapped_allalign.bam
    """
}
