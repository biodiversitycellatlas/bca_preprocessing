process FILTER_MATRICES_ALEVIN {
    publishDir "${params.outdir}/mapping_alevin/${meta.id}/${meta.id}_counts/alevin/filtered_secondderiv", mode: 'copy'
    label 'process_single_mem2'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(alevin_mtx_dir), path(cutoff_txt)

    output:
    tuple val(meta), path("filtered/"),                                 emit: filtered_matrix
    tuple val(meta), path("${meta.id}_secondderiv_statistics.json"),    emit: filter_stats

    script:
    def usa_counts = params.alevin_usa_counts ?: 'SUA'
    """
    # Read the integer cutoff from file
    CUTOFF=\$(cat ${cutoff_txt})

    echo "Applying UMI threshold of \$CUTOFF to ${meta.id}"

    secondderiv_alevin.py filter \\
        -d ${alevin_mtx_dir} \\
        -c \$CUTOFF \\
        -o filtered/ \\
        -s ${meta.id}_secondderiv_statistics.json \\
        --counts ${usa_counts}
    """
}
