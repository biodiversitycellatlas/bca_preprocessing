process CALC_MT_RRNA {
    publishDir "${params.outdir}/rRNA_mtDNA", mode: 'copy'
    label 'process_single'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.1.1--h577a1d6_0' :
        'biocontainers/subread:2.1.1--h577a1d6_0' }"

    input:
    tuple val(meta), path(mapping_files)
    file(bam_index)

    output:
    file("${meta.id}_mt_rrna_metrics.txt")

    script:
    """
    echo "\n\n==================  CALCULATION rRNA & mtDNA =================="
    echo "Sample ID: ${meta}"

    bam_file=\$(ls ${meta.id}_Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"

    calculate_rrna_mtdna.sh \\
        \${bam_file} \\
        ${meta.id}_mt_rrna_metrics.txt \\
        ${params.ref_gtf} \\
        ${params.grep_rrna} \\
        ${params.mt_contig}

    echo "Created file: ${meta.id}_mt_rrna_metrics.txt"
    """
}
