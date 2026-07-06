process FILTER_MATRICES {
    publishDir "${params.outdir}/mapping_STARsolo/${meta.id}/${meta.id}_Solo.out/GeneFull_Ex50pAS/filtered_secondderiv", mode: 'copy'
    label 'process_low'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(raw_matrix_dir)
    path(cutoff_txt)

    output:
    tuple val(meta), path("filtered/"),                     emit: filtered_matrix
    tuple val(meta), path("secondderiv_statistics.json"),   emit: filter_stats

    script:
    """
    # Read the integer cutoff from file
    CUTOFF=\$(cat ${cutoff_txt})

    echo "Applying UMI threshold of \$CUTOFF to ${meta.id}"

    secondderiv_filter_matrices.py \\
        --dir ${raw_matrix_dir} \\
        --cutoff \$CUTOFF \\
        --outdir filtered/
    """
}
