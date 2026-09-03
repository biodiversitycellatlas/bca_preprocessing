process DOUBLET_FILTER {
    publishDir "${params.outdir}/doublet_filtering/${meta.id}/${meta.datatype}", mode: 'copy'
    tag "${meta.id} | ${meta.mapping_method} | ${meta.datatype}"
    label 'process_low'

    // Memory tracks the size of the matrix being sliced, overrides process_low's flat
    // assignments. Coefficients live in params.dynamic_memory; remove the entry to fall
    // back to the plain label.
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.DOUBLET_FILTER, [mtx, barcodes, features], task.attempt, 12) }

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(mtx), path(barcodes), path(features), path(combined_results)

    output:
    tuple val(meta), path("${meta.id}_${meta.datatype}_doublet_filtered"),      emit: matrix
    path("${meta.id}_${meta.datatype}_doublet_filter_summary.txt"),             emit: summary
    path("${meta.id}_${meta.datatype}_doublet_filtering_summary.png"),          emit: filtering_summary_plot

    script:
    """
    echo "\n\n==================  Doublet consensus filter =================="
    echo "Meta: ${meta}"
    echo "Matrix: ${mtx}"
    echo "Combined results: ${combined_results}"

    filter_doublets.py \\
        --mtx ${mtx} \\
        --barcodes ${barcodes} \\
        --features ${features} \\
        --combined_results ${combined_results} \\
        --method ${params.doublet_consensus_method} \\
        --outdir ${meta.id}_${meta.datatype}_doublet_filtered \\
        --summary_txt ${meta.id}_${meta.datatype}_doublet_filter_summary.txt \\
        --image_prefix ${meta.id}_${meta.datatype}_
    """
}
