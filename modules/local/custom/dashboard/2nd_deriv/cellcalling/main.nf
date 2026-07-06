process SECONDDERIV_CELLCALLING {
    label 'process_low'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(umiPerCell_sorted)

    output:
    tuple val(meta), path("${meta.id}_knee_data.json"),     emit: json_data
    path("${meta.id}_cutoff.txt"),                          emit: cutoff

    script:
    """
    secondderiv_cellcalling.py \\
        -i ${umiPerCell_sorted} \\
        -s ${meta.id} \\
        -o ${meta.id}_knee_data.json \\
        -c ${meta.id}_cutoff.txt
    """
}
