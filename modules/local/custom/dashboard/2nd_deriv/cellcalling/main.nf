process SECONDDERIV_CELLCALLING {
    publishDir path: { mapper == 'alevin-fry'
                        ? "${params.outdir}/mapping_alevin/${meta.id}/${meta.id}_counts/alevin/filtered_secondderiv"
                        : "${params.outdir}/mapping_STARsolo/${meta.id}/${meta.id}_Solo.out/GeneFull_Ex50pAS/filtered_secondderiv" },
               mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${mapper}"

    // STARsolo hands over its own UMI curve, a small text file (process_single2's 500 MB);
    // alevin-fry's curve is derived here from the whole quantified matrix (process_single_mem2's
    // 10 GB). Set here rather than through a label: a withLabel block in the config takes
    // precedence over a directive in the process, so a label would pin both mappers to one size.
    cpus   1
    memory { (mapper == 'alevin-fry' ? 10.GB : 500.MB) * task.attempt }
    time   { 30.m * task.attempt }

    conda "${moduleDir}/environment.yml"

    input:
    // mapper: 'starsolo' or 'alevin-fry'. curve_source: STARsolo's UMIperCellSorted.txt, or
    // alevin-fry's <sample>_counts/alevin directory, from which the curve is derived
    tuple val(meta), val(mapper), path(curve_source)

    output:
    tuple val(meta), path("${meta.id}_alevin_UMIperCellSorted.txt"), emit: umi_per_cell, optional: true
    tuple val(meta), path("${meta.id}_knee_data.json"),              emit: json_data
    tuple val(meta), path("${meta.id}_cutoff.txt"),                  emit: cutoff
    path "versions.yml",                                             emit: versions

    script:
    def usa_counts = params.alevin_usa_counts ?: 'SUA'
    // With "manual_cutoff" the threshold is taken from the samplesheet and the search is
    // skipped, so expected_cells (which only bounds that search) has nothing left to do.
    // The UMI curve is still needed, since the report plots it either way.
    def manual_cutoff = params.cellfilter_method == "manual_cutoff" ? meta.manual_cutoff : null
    if (params.cellfilter_method == "manual_cutoff" && manual_cutoff == null) {
        error "Sample '${meta.id}' has no 'manual_cutoff' in the samplesheet, which cellfilter_method = 'manual_cutoff' requires"
    }
    def cutoff_arg = manual_cutoff != null ? "-m ${manual_cutoff}" : ""
    def expected_cells_arg = (manual_cutoff == null && meta.expected_cells) ? "-e ${meta.expected_cells}" : ""

    // alevin-fry has no UMI-per-cell curve of its own, so one is written in STARsolo's format first
    def umi_curve = mapper == 'alevin-fry' ? "${meta.id}_alevin_UMIperCellSorted.txt" : "${curve_source}"
    def derive_curve = mapper == 'alevin-fry'
        ? "secondderiv_alevin.py umis -d ${curve_source} -o ${umi_curve} --counts ${usa_counts}"
        : ""
    """
    ${derive_curve}

    secondderiv_cellcalling.py \\
        -i ${umi_curve} \\
        -s ${meta.id} \\
        -o ${meta.id}_knee_data.json \\
        -c ${meta.id}_cutoff.txt \\
        ${expected_cells_arg} ${cutoff_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
