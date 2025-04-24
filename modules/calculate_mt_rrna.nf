process CALC_MT_RRNA {   
    publishDir "${params.output_dir}/rRNA_mtDNA", mode: 'copy'
    debug true

    input:
    tuple val(meta), path(mapping_files)
    file(bam_index)

    output:
    file("${meta.id}_mt_rrna_metrics.txt")     
       
    script:
    """
    echo "\n\n==================  CALCULATION rRNA & mtDNA =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Sample ID: ${meta}"

    bam_file=\$(ls ${meta.id}_Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"

    bash ${launchDir}/bin/calculate_rrna_mtdna.sh \\
        \${bam_file} \\
        ${meta.id}_mt_rrna_metrics.txt \\
        ${params.ref_gtf} \\
        ${params.grep_rrna} \\
        ${params.mt_contig}

    echo "Created file: ${meta.id}_mt_rrna_metrics.txt"
    """
}
