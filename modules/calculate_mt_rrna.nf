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
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Sample ID: ${sample_id}"
    echo "Config: ${config}"

    bam_file=\$(ls ${sample_id}_Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"

    sbatch ${params.baseDir}/bin/calculate_rrna_mtdna.sh \\
        \${bam_file} \\
        ${sample_id}_${config}_mt_rrna_metrics.txt \\
        ${params.ref_star_gtf} \\
        ${params.grep_rrna} \\
        ${params.mt_contig}
    """
}
