process FILTER_MATRICES {
    publishDir "${params.outdir}/mapping_STARsolo/${meta.id}/${meta.id}_Solo.out/GeneFull_Ex50pAS/filtered_secondderiv", mode: 'copy'
    label 'process_single_mem2'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(raw_matrix_dir), path(cutoff_txt), path(cellreads_stats)

    output:
    tuple val(meta), path("filtered/"),                                 emit: filtered_matrix
    tuple val(meta), path("${meta.id}_secondderiv_statistics.json"),    emit: filter_stats

    script:
    // CellReads.stats is optional: without it the read-level statistics are skipped
    def cellreads_arg = cellreads_stats ? "--cellreads ${cellreads_stats}" : ""
    """
    # Read the integer cutoff from file
    CUTOFF=\$(cat ${cutoff_txt})

    echo "Applying UMI threshold of \$CUTOFF to ${meta.id}"

    secondderiv_filter_matrices.py \\
        --dir ${raw_matrix_dir} \\
        --cutoff \$CUTOFF \\
        --outdir filtered/ \\
        --stats ${meta.id}_secondderiv_statistics.json \\
        ${cellreads_arg}
    """
}
