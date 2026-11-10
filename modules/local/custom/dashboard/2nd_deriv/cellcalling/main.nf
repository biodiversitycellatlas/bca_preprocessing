process SECONDDERIV_CELLCALLING {
    publishDir path: { mapper == 'alevin-fry'
                        ? "${params.outdir}/mapping_alevin/${meta.id}/${meta.id}_counts/alevin/filtered_secondderiv"
                        : "${params.outdir}/mapping_STARsolo/${meta.id}/${meta.id}_Solo.out/GeneFull_Ex50pAS/filtered_secondderiv" },
               mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${mapper}"

    // The caller reads the raw count matrix (its ambient test needs it), and the curve is
    // derived here from that same matrix, so memory tracks the matrix (only the analysed one
    // in STARsolo's raw/ is measured).
    cpus   1
    memory {
        BcaResources.scaledMemory(
            params.dynamic_memory?.SECONDDERIV_CELLCALLING,
            mapper == 'alevin-fry'
                ? [curve_source]
                : BcaResources.filesIn(curve_source, [matrix_name, 'barcodes.tsv', 'features.tsv']),
            task.attempt, 10)
    }
    time   { 30.m * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_pandas_scipy:eeb34b63825ecb6a' :
        'community.wave.seqera.io/library/python_numpy_pandas_scipy:b12b357472dfbd06' }"

    input:
    // mapper: 'starsolo' or 'alevin-fry'. curve_source: the raw count directory, STARsolo's
    // GeneFull_Ex50pAS/raw or alevin-fry's <sample>_counts/alevin. matrix_name: the analysed
    // matrix in STARsolo's raw/ (StarsoloMatrix.rawFileName), '' for alevin-fry.
    tuple val(meta), val(mapper), path(curve_source), val(matrix_name)

    output:
    // <id>_alevin_UMIperCellSorted.txt, or for STARsolo <id>_<multi-gene type>_ or <id>_Unique_
    tuple val(meta), path("${meta.id}_*_UMIperCellSorted.txt"),     emit: umi_per_cell, optional: true
    tuple val(meta), path("${meta.id}_knee_data.json"),              emit: json_data
    tuple val(meta), path("${meta.id}_cutoff.txt"),                  emit: cutoff
    path "versions.yml",                                             emit: versions

    script:
    def usa_counts = params.alevin_usa_counts ?: 'SUA'
    
    // With "manual_cutoff" the threshold is taken from the samplesheet and the search is
    // skipped, so expected_cells (which only informs that search) has nothing left to do.
    // The UMI curve is still needed, since the report plots it either way.
    def manual_cutoff = params.cellfilter_method == "manual_cutoff" ? meta.manual_cutoff : null
    if (params.cellfilter_method == "manual_cutoff" && manual_cutoff == null) {
        error "Sample '${meta.id}' has no 'manual_cutoff' in the samplesheet, which cellfilter_method = 'manual_cutoff' requires"
    }
    def cutoff_arg = manual_cutoff != null ? "-m ${manual_cutoff}" : ""
    def expected_cells_arg = (manual_cutoff == null && meta.expected_cells) ? "-e ${meta.expected_cells}" : ""

    // The curve is written in STARsolo's format from the analysed matrix, so that it and the
    // matrix the caller tests give the same barcode ranks
    def multi_match = matrix_name ? (matrix_name =~ /^UniqueAndMult-(.+)\.mtx$/) : null
    def multi_type = multi_match ? multi_match[0][1] : 'Unique'
    def umi_curve = mapper == 'alevin-fry' ? "${meta.id}_alevin_UMIperCellSorted.txt"
        : "${meta.id}_${multi_type}_UMIperCellSorted.txt"
    def derive_curve = mapper == 'alevin-fry'
        ? "secondderiv_alevin.py umis -d ${curve_source} -o ${umi_curve} --counts ${usa_counts}"
        : "secondderiv_filter_matrices.py umis -d ${curve_source} -m ${matrix_name} -o ${umi_curve}"
    def matrix_arg = mapper == 'alevin-fry' ? "--counts ${usa_counts}" : "--matrix ${matrix_name}"
    """
    ${derive_curve}

    secondderiv_cellcalling.py \\
        -i ${umi_curve} \\
        -s ${meta.id} \\
        -o ${meta.id}_knee_data.json \\
        -c ${meta.id}_cutoff.txt \\
        --raw-dir ${curve_source} ${matrix_arg} \\
        ${expected_cells_arg} ${cutoff_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
