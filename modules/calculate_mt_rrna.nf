// =============  CALCULATION rRNA & mtRNA  ============= \\ 
// 

process CALC_MT_RRNA {   
    publishDir "${params.resDir}/rRNA_mtDNA/rRNA_mtDNA_${config_name}/${sample_id}", mode: 'copy'
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
    samtools view -h -F 256 Aligned.sortedByCoord.out.bam | grep -E "^\\@|NH:i:[2-9]" | samtools view -b -o multimapped_primealign.bam
    # Create seperate BAM with filtered results - filter for multimapped reads present in both secondary and primary alignments
    samtools view -h Aligned.sortedByCoord.out.bam | grep -E "^\\@|NH:i:[2-9]" | samtools view -b -o multimapped_allalign.bam

    # Calculation rRNA
    featureCounts -t rRNA -a ${params.ref_star_gtf} -o feat_counts_rRNA.txt Aligned.sortedByCoord.out.bam
    featureCounts -M --fraction -t rRNA -a ${params.ref_star_gtf} -o feat_counts_rRNA_mmpa.txt multimapped_primealign.bam
    featureCounts -M --fraction -t rRNA -a ${params.ref_star_gtf} -o feat_counts_rRNA_mmaa.txt multimapped_allalign.bam

    # Calculation mtDNA
    grep 'mtDNA' ${params.ref_star_gtf} > mtDNA_only.gtf
    featureCounts -t exon -a mtDNA_only.gtf -o feat_counts_mtDNA.txt Aligned.sortedByCoord.out.bam
    featureCounts -M --fraction -t exon -a mtDNA_only.gtf -o feat_counts_mtDNA_mmpa.txt multimapped_primealign.bam
    featureCounts -M --fraction -t exon -a mtDNA_only.gtf -o feat_counts_mtDNA_mmaa.txt multimapped_allalign.bam

    # Per-gene counts for multimappers
    featureCounts -M --fraction -a ${params.ref_star_gtf} -o feat_counts_multimappers_genelevel.txt multimapped_allalign.bam
    
    # Saves feature counts as clean gene counts table
    tail -n +3 feat_counts_multimappers_genelevel.txt | cut -f1,7 > multimapper_gene_counts.txt
    
    # Sort the counts
    sort -k2,2nr multimapper_gene_counts.txt > sorted_multimapper_gene_counts.txt
    """
}
