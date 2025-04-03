process CALC_MT_RRNA {   
    publishDir "${params.resDir}/rRNA_mtDNA/${config}", mode: 'copy'
    debug true

    input:
    tuple val(sample_id), path(mapping_files)
    file(bam_index)
    val(config)

    output:
    path("*")     
       
    script:
    """
    echo "\n\n==================  CALCULATION rRNA & mtDNA =================="
    echo "Sample ID: ${sample_id}"
    echo "Config: ${config}"

    sbatch ${params.baseDir}/bin/calculate_rrna_mtdna.sh \\
        ${params.resDir} \\
        ${sample_id} \\
        ${config} \\
        ${params.ref_star_gtf} \\
        ${params.grep_rrna} \\
        ${params.mt_contig}
    """
}
