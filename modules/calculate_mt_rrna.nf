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

    # Calculation rRNA
    featureCounts \\
        -t rRNA \\
        -a ${params.ref_star_gtf} \\
        -o feat_counts_rRNA.txt \\
        Aligned.sortedByCoord.out.bam

    # Calculation mtDNA
    grep 'mtDNA' ${params.ref_star_gtf} > mtDNA_only.gtf
    featureCounts \\
        -t exon \\
        -a mtDNA_only.gtf \\
        -o feat_counts_mtDNA.txt \\
        Aligned.sortedByCoord.out.bam
    """
}
