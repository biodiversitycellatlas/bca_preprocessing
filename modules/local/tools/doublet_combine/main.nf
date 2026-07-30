process COMBINE_DOUBLET_RESULTS {
    publishDir "${params.outdir}/doublet_filtering/${meta.id}/combined", mode: 'copy'
    tag "${meta.id} | ${mapping_method}"
    label 'process_single'

    container { demuxafy_sif }

    input:
    tuple val(meta), path(scrublet_results), path(scdblfinder_results), val(mapping_method), val(datatype)
    val demuxafy_sif

    output:
    tuple val(meta), path("${meta.id}_combined_doublets_w_combined_assignments.tsv"), val(mapping_method), val(datatype), emit: combined_results
    path("${meta.id}_combined_doublets_summary.tsv"),                                                                    emit: combined_summary, optional: true

    script:
    """
    echo "\n\n==================  Combine doublet results =================="
    echo "Meta: ${meta}"
    echo "Consensus method: ${params.doublet_consensus_method}"

    Combine_Results.R \\
        -o ${meta.id}_combined_doublets.tsv \\
        -r ${scrublet_results} \\
        -n ${scdblfinder_results} \\
        -m "${params.doublet_consensus_method}"
    """
}
