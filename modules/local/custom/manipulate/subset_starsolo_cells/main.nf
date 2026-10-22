process SUBSET_STARSOLO_CELLS {
    // Published as e.g. GeneFull_Ex50pAS/filtered_UniqueAndMult-EM, next to STAR's own filtered/ (unique counts)
    publishDir path: { "${params.outdir}/mapping_STARsolo/${meta.id}/${meta.id}_Solo.out/GeneFull_Ex50pAS" }, mode: 'copy',
               saveAs: { filename -> filename == 'versions.yml' ? null
                                   : filename == 'filtered' ? "filtered_${matrix_name.replaceAll(/\.mtx$/, '')}"
                                   : filename }
    tag "${meta.id}"
    label 'process_single_mem2'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_pandas_scipy:eeb34b63825ecb6a' :
        'community.wave.seqera.io/library/python_numpy_pandas_scipy:b12b357472dfbd06' }"

    input:
    // raw_dir: STARsolo's GeneFull_Ex50pAS/raw. cell_barcodes: STAR's own call (filtered/barcodes.tsv).
    // matrix_name: the multi-gene matrix in raw_dir to subset (StarsoloMatrix.rawFileName)
    tuple val(meta), path(raw_dir), path(cell_barcodes), val(matrix_name)

    output:
    tuple val(meta), path("filtered/"),                    emit: filtered_matrix
    tuple val(meta), path("${meta.id}_*_statistics.json"), emit: matrix_stats
    path "versions.yml",                                   emit: versions

    script:
    // e.g. <id>_UniqueAndMult-EM_statistics.json, published next to filtered_UniqueAndMult-EM/
    def stats_json = "${meta.id}_${matrix_name.replaceAll(/\.mtx$/, '')}_statistics.json"
    """
    echo "Subsetting ${matrix_name} of ${meta.id} to STARsolo's \$(wc -l < ${cell_barcodes}) called cells"

    subset_matrices_to_cells.py \\
        --dir ${raw_dir} \\
        --barcodes ${cell_barcodes} \\
        --matrices ${matrix_name} \\
        --outdir filtered/

    # Every matrix directory the pipeline writes holds the analysed counts as matrix.mtx
    mv filtered/${matrix_name} filtered/matrix.mtx

    # STAR's Summary.csv describes the unique-gene counts of these cells, so the dashboard
    # takes the median UMIs, median genes and genes detected of this matrix from here
    secondderiv_filter_matrices.py stats \\
        --dir filtered/ \\
        --matrix matrix.mtx \\
        --label ${matrix_name} \\
        --stats ${stats_json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
