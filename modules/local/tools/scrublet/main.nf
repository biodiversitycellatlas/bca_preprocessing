process SCRUBLET {
    publishDir path: { "${params.outdir}/doublet_filtering/${meta.id}/scrublet" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${meta.mapping_method}"

    // Memory tracks the cell-called matrix
    cpus   1
    time   { BcaResources.attemptTime(task, '45m') }
    memory { BcaResources.attemptMemory(task, BcaResources.scaledMemory(
        params.dynamic_memory?.SCRUBLET, [tenx_dir], 1, 6)) }

    // Scrublet cannot model every matrix (too few cells, no bimodal score split). A
    // failure never ends the run (error_optional: one retry when the task was killed for
    // its resources, otherwise ignored); the sample's doublets then come from scDblFinder
    // alone, and filtering_workflow.nf says so.
    label 'error_optional'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_scrublet_matplotlib-base_numpy_pruned:92c538727b8dc48e' :
        'community.wave.seqera.io/library/python_scrublet_matplotlib-base_numpy_pruned:22ac1f4745d87ed0' }"

    input:
    tuple val(meta), path(tenx_dir)

    output:
    tuple val(meta), path("${meta.id}_scrublet_results.tsv"), emit: scrublet_results
    path("${meta.id}_doublet_score_histogram.png"),           emit: scrublet_histogram, optional: true
    path "versions.yml",                                      emit: versions

    script:
    """
    echo -e "\\n\\n==================  Scrublet =================="
    echo "Meta: ${meta}"
    echo "10x dir: ${tenx_dir}"
    echo "Datatype: ${meta.datatype}"

    # numba (pulled in by scrublet) and matplotlib cache to a writable directory, not
    # into the read-only image
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    export MPLCONFIGDIR=\$PWD/.mplconfig

    run_scrublet_doublet.py \\
        --counts_matrix ${tenx_dir} \\
        --outdir scrublet_out \\
        --sample_id ${meta.id} \\
        --doublet_rate_per_1000 ${params.doublet_rate_per_1000}

    mv scrublet_out/scrublet_results.tsv ${meta.id}_scrublet_results.tsv
    for histogram in scrublet_out/doublet_score_histogram_manual_threshold.png scrublet_out/doublet_score_histogram.png; do
        if [ -f "\$histogram" ]; then
            mv "\$histogram" ${meta.id}_doublet_score_histogram.png
            break
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        scrublet: \$(python3 -c "import importlib.metadata as m; print(m.version('scrublet'))")
    END_VERSIONS
    """
}
