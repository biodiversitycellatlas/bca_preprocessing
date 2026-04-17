process CELLSWEEP {
    publishDir "${params.outdir}/cellsweep/${meta.id}", mode: 'copy'
    tag "${meta.id} | ${mapping_method}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // container "oras://community.wave.seqera.io/library/anndata_numpy_pandas_python_pruned:a3b0b95a49665473"

    input:
    tuple val(meta), path(raw_h5ad), val(mapping_method), val(datatype)

    output:
    path("${meta.id}_cs_filtered.h5ad"),      emit: cs_filtered_h5ad
    path("${meta.id}_cs_full.h5ad"),          emit: cs_full_h5ad
    path("*_ambient_hat_histogram.png"),      emit: cs_ambient_hist_plot
    path("*_top_ambient_genes.csv"),          emit: cs_top_genes
    path("*_umap_comparison.png"),            emit: cs_umap_comparison_plot
    path("*_doublet_summary.txt"),            emit: cs_doublet_summary, optional: true

    script:
    """
    echo "\n\n==================  CellSweep =================="
    echo "Meta: ${meta}"
    echo "Raw h5ad: ${raw_h5ad}"
    echo "Mapping method: ${mapping_method}"
    echo "Datatype: ${datatype}"

    run_cellsweep.py \\
        --input_h5ad ${raw_h5ad} \\
        --cs_filtered_h5ad ${meta.id}_cs_filtered.h5ad \\
        --cs_full_h5ad ${meta.id}_cs_full.h5ad \\
        --image_prefix ${meta.id}_ \\
        --expected_cells ${meta.expected_cells} \\
        --threads ${task.cpus}
    """
}
