process SCRUBLET {
    publishDir "${params.outdir}/doublet_filtering/${meta.id}/scrublet", mode: 'copy'
    tag "${meta.id} | ${mapping_method}"
    label 'process_medium'

    container { demuxafy_sif }

    input:
    tuple val(meta), path(tenx_dir), val(mapping_method), val(datatype), path(filtered_barcodes)
    val demuxafy_sif

    output:
    tuple val(meta), path("${meta.id}_scrublet_results.tsv"), val(mapping_method), val(datatype), emit: scrublet_results
    path("${meta.id}_doublet_score_histogram.png"),                                              emit: scrublet_histogram, optional: true
    path "versions.yml",                                                                          emit: versions

    script:
    // STARsolo's raw output includes near-empty droplets that break Scrublet's automatic
    // doublet-rate estimate; restrict doublet calling to real cells when a filtered list is given.
    def filter_arg = filtered_barcodes.name != 'NO_FILE' ? "--filtered_barcodes ${filtered_barcodes}" : ""
    """
    echo "\n\n==================  Scrublet =================="
    echo "Meta: ${meta}"
    echo "10x dir: ${tenx_dir}"
    echo "Filtered barcodes: ${filtered_barcodes}"

    run_scrublet_doublet.py \\
        --counts_matrix ${tenx_dir} \\
        --outdir scrublet_out \\
        --sample_id ${meta.id} \\
        ${filter_arg}

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
