process SECONDDERIV_CELLCALLING_ALEVIN {
    publishDir "${params.outdir}/mapping_alevin/${meta.id}/${meta.id}_counts/alevin/filtered_secondderiv", mode: 'copy'
    label 'process_single2'
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
    // With "manual_cutoff" the threshold is taken from the samplesheet and the search is
    // skipped, so expected_cells (which only bounds that search) has nothing left to do.
    // The UMI curve is still derived, since the report plots it either way.
    def manual_cutoff = params.cellfilter_method == "manual_cutoff" ? meta.manual_cutoff : null
    if (params.cellfilter_method == "manual_cutoff" && manual_cutoff == null) {
        error "Sample '${meta.id}' has no 'manual_cutoff' in the samplesheet, which cellfilter_method = 'manual_cutoff' requires"
    }
    def cutoff_arg = manual_cutoff != null ? "-m ${manual_cutoff}" : ""
    def expected_cells_arg = (manual_cutoff == null && meta.expected_cells) ? "-e ${meta.expected_cells}" : ""
    """
    secondderiv_alevin.py umis \\
        -d ${alevin_mtx_dir} \\
        -o ${meta.id}_alevin_UMIperCellSorted.txt \\
        --counts ${usa_counts}

    secondderiv_cellcalling.py \\
        -i ${meta.id}_alevin_UMIperCellSorted.txt \\
        -s ${meta.id} \\
        -o ${meta.id}_knee_data.json \\
        -c ${meta.id}_cutoff.txt \\
        ${expected_cells_arg} ${cutoff_arg}
    """
}
