process CELLSWEEP {
    publishDir "${params.outdir}/cellsweep/${meta.id}", mode: 'copy'
    tag "${meta.id} | ${meta.mapping_method}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // container "oras://community.wave.seqera.io/library/anndata_numpy_pandas_python_pruned:a3b0b95a49665473"

    input:
    tuple val(meta), path(raw_h5ad), path(doublet_results)

    output:
    path("${meta.id}_cs_filtered.h5ad"),      emit: cs_filtered_h5ad
    path("${meta.id}_cs_full.h5ad"),          emit: cs_full_h5ad
    path("*_ambient_hat_histogram.png"),      emit: cs_ambient_hist_plot
    path("*_top_ambient_genes.csv"),          emit: cs_top_genes
    path("*_umap_comparison.png"),            emit: cs_umap_comparison_plot
    path "versions.yml",                      emit: versions

    script:
    def doublet_arg = doublet_results ? "--doublet_results ${doublet_results} --doublet_method ${params.doublet_consensus_method}" : ""
    """
    echo "\n\n==================  CellSweep =================="
    echo "Meta: ${meta}"
    echo "Raw h5ad: ${raw_h5ad}"
    echo "Mapping method: ${meta.mapping_method}"
    echo "Datatype: ${meta.datatype}"
    echo "Doublet results: ${doublet_results ?: 'none'}"

    run_cellsweep.py \\
        --input_h5ad ${raw_h5ad} \\
        --cs_filtered_h5ad ${meta.id}_cs_filtered.h5ad \\
        --cs_full_h5ad ${meta.id}_cs_full.h5ad \\
        --image_prefix ${meta.id}_ \\
        --expected_cells ${meta.expected_cells} \\
        --threads ${task.cpus} \\
        ${doublet_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
