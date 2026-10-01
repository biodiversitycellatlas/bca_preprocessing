process CELLSWEEP {
    publishDir "${params.outdir}/cellsweep/${meta.id}", mode: 'copy'
    tag "${meta.id} | ${meta.mapping_method}"
    label 'process_low2'
    debug true

    // A matrix with too few cells to cluster exits 3 from bin/run_cellsweep.py; the
    // sample is then dropped from this step and published without denoised counts.
    // Other failures keep the error_optional rule from conf/base.config: one retry for
    // the transient exit codes, otherwise skipped the same way.
    errorStrategy {
        if (task.exitStatus in ((130..145) + 104 + 175 + 255)) {
            return 'retry'
        }
        def reason = task.exitStatus == 3
            ? "too few cells passed its filters"
            : "it failed with exit status ${task.exitStatus}"
        log.warn "No CellSweep output for ${meta.id} (${meta.mapping_method}): ${reason} " +
            "(see .command.log in the task work directory). The sample continues without denoised counts."
        return 'ignore'
    }
    maxRetries 1

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(mtx), path(barcodes), path(features), path(doublet_results)

    output:
    tuple val(meta), path("${meta.id}_cs_filtered.h5ad"), emit: cs_filtered_h5ad
    tuple val(meta), path("${meta.id}_cs_full.h5ad"),     emit: cs_full_h5ad
    path("*_ambient_hat_histogram.png"),                  emit: cs_ambient_hist_plot
    path("*_top_ambient_genes.csv"),                      emit: cs_top_genes
    path("*_umap_comparison.png"),                        emit: cs_umap_comparison_plot
    path "versions.yml",                                  emit: versions

    script:
    def doublet_arg = doublet_results ? "--doublet_results ${doublet_results} --doublet_method ${params.doublet_consensus_method}" : ""
    """
    echo "\n\n==================  CellSweep =================="
    echo "Meta: ${meta}"
    echo "Raw matrix: ${mtx}"
    echo "Mapping method: ${meta.mapping_method}"
    echo "Datatype: ${meta.datatype}"
    echo "Doublet results: ${doublet_results ?: 'none (already filtered, or detection off)'}"

    run_cellsweep.py \\
        --mtx ${mtx} \\
        --barcodes ${barcodes} \\
        --features ${features} \\
        --sample_id ${meta.id} \\
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
