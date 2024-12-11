// =============  CALCULATION rRNA & mtRNA  ============= \\ 
// 

process CALC_MT_RRNA {   
    publishDir "{params.resDir}/mapping_STARsolo_${config_name}/${sample_id}_rRNA", mode: 'copy'
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    output:
    path("*")     
       
    script:
    """
    echo "\n\n==================  CALCULATION rRNA & mtRNA ${config_name} =================="
    echo "Sample ID: ${sample_id}"
    echo "Config name: ${config_name}"
    echo "Mapping files: ${mapping_files}"

    # Calculation rRNA
    featureCounts -M --fraction -f \\
        -t rRNA \\
        -a ${params.gtf_file} \\
        -o feat_counts_rRNA.txt \\
        ${mapping_files}/*Aligned.sortedByCoord.out.bam
    """
}
