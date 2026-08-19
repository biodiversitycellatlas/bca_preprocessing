process SCRUBLET {
    publishDir "${params.outdir}/doublet_filtering/${meta.id}/scrublet", mode: 'copy'
    tag "${meta.id} | ${meta.mapping_method}"
    label 'process_low'

    // Scrublet cannot model every matrix (too few cells, no bimodal score split). The
    // sample then continues without doublet annotation instead of ending the run.
    label 'error_optional'

    container { demuxafy_sif }

    input:
    tuple val(meta), path(tenx_dir)
    val demuxafy_sif

    output:
    tuple val(meta), path("${meta.id}_scrublet_results.tsv"), emit: scrublet_results
    path("${meta.id}_doublet_score_histogram.png"),           emit: scrublet_histogram, optional: true
    path "versions.yml",                                      emit: versions

    script:
    """
    echo "\n\n==================  Scrublet =================="
    echo "Meta: ${meta}"
    echo "10x dir: ${tenx_dir}"
    echo "Datatype: ${meta.datatype}"

    run_scrublet_doublet.py \\
        --counts_matrix ${tenx_dir} \\
        --outdir scrublet_out \\
        --sample_id ${meta.id}

    mv scrublet_out/scrublet_results.tsv ${meta.id}_scrublet_results.tsv
    if [ -f scrublet_out/doublet_score_histogram.png ]; then
        mv scrublet_out/doublet_score_histogram.png ${meta.id}_doublet_score_histogram.png
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
