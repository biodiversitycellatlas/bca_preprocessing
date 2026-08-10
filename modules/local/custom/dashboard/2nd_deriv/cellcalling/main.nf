process SECONDDERIV_CELLCALLING {
    publishDir "${params.outdir}/mapping_STARsolo/${meta.id}/${meta.id}_Solo.out/GeneFull_Ex50pAS/filtered_secondderiv", mode: 'copy'
    label 'process_low'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(umiPerCell_sorted)

    output:
    tuple val(meta), path("${meta.id}_knee_data.json"),     emit: json_data
    tuple val(meta), path("${meta.id}_cutoff.txt"),         emit: cutoff

    script:
    def expected_cells_arg = meta.expected_cells ? "-e ${meta.expected_cells}" : ""
    """
    secondderiv_cellcalling.py \\
        -i ${umiPerCell_sorted} \\
        -s ${meta.id} \\
        -o ${meta.id}_knee_data.json \\
        -c ${meta.id}_cutoff.txt \\
        ${expected_cells_arg}
    """
}
