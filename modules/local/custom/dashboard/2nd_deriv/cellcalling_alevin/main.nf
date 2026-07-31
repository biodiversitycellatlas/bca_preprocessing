process SECONDDERIV_CELLCALLING_ALEVIN {
    publishDir "${params.outdir}/mapping_alevin/${meta.id}/${meta.id}_counts/alevin/filtered_secondderiv", mode: 'copy'
    label 'process_low'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(alevin_mtx_dir)

    output:
    tuple val(meta), path("${meta.id}_alevin_UMIperCellSorted.txt"), emit: umi_per_cell
    tuple val(meta), path("${meta.id}_knee_data.json"),              emit: json_data
    tuple val(meta), path("${meta.id}_cutoff.txt"),                  emit: cutoff

    script:
    def usa_counts = params.alevin_usa_counts ?: 'SUA'
    """
    secondderiv_alevin.py umis \\
        -d ${alevin_mtx_dir} \\
        -o ${meta.id}_alevin_UMIperCellSorted.txt \\
        --counts ${usa_counts}

    secondderiv_cellcalling.py \\
        -i ${meta.id}_alevin_UMIperCellSorted.txt \\
        -s ${meta.id} \\
        -o ${meta.id}_knee_data.json \\
        -c ${meta.id}_cutoff.txt
    """
}
