process FILTER_MATRICES {
    publishDir path: { mapper == 'alevin-fry'
                        ? "${params.outdir}/mapping_alevin/${meta.id}/${meta.id}_counts/alevin/filtered_secondderiv"
                        : "${params.outdir}/mapping_STARsolo/${meta.id}/${meta.id}_Solo.out/GeneFull_Ex50pAS/filtered_secondderiv" },
               mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_single_mem2'
    tag "${meta.id} | ${mapper}"

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/python_numpy_pandas_scipy:eeb34b63825ecb6a"

    input:
    // mapper: 'starsolo' or 'alevin-fry'. matrix_dir: STARsolo's GeneFull_Ex50pAS/raw, or
    // alevin-fry's <sample>_counts/alevin. cellreads_stats: STARsolo's CellReads.stats, [] otherwise
    tuple val(meta), val(mapper), path(matrix_dir), path(cutoff_txt), path(cellreads_stats)

    output:
    tuple val(meta), path("filtered/"),                                 emit: filtered_matrix
    tuple val(meta), path("${meta.id}_secondderiv_statistics.json"),    emit: filter_stats
    path "versions.yml",                                                emit: versions

    script:
    // The two mappers lay their matrices out differently (genes x cells vs cells x genes, with
    // alevin-fry's USA blocks side by side), so each has its own filter script.
    // CellReads.stats is optional: without it the read-level statistics are skipped
    def cellreads_arg = cellreads_stats ? "--cellreads ${cellreads_stats}" : ""
    def usa_counts = params.alevin_usa_counts ?: 'SUA'
    def filter_cmd = mapper == 'alevin-fry'
        ? "secondderiv_alevin.py filter -d ${matrix_dir} -c \$CUTOFF -o filtered/ -s ${meta.id}_secondderiv_statistics.json --counts ${usa_counts}"
        : "secondderiv_filter_matrices.py --dir ${matrix_dir} --cutoff \$CUTOFF --outdir filtered/ --stats ${meta.id}_secondderiv_statistics.json ${cellreads_arg}"
    """
    # Read the integer cutoff from file
    CUTOFF=\$(cat ${cutoff_txt})

    echo "Applying UMI threshold of \$CUTOFF to ${meta.id}"

    ${filter_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """
}
