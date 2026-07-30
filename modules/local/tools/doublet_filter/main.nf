process DOUBLET_FILTER {
    publishDir "${params.outdir}/doublet_filtering/${meta.id}/filtered", mode: 'copy'
    tag "${meta.id} | ${mapping_method}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(raw_h5ad), path(combined_results), val(mapping_method), val(datatype)

    output:
    tuple val(meta), path("${meta.id}_doublet_filtered.h5ad"), val(mapping_method), val(datatype), emit: h5ad
    path("${meta.id}_doublet_filter_summary.txt"),                                                emit: summary
    path("${meta.id}_doublet_filtering_summary.png"),                                             emit: filtering_summary_plot

    script:
    """
    echo "\n\n==================  Doublet consensus filter =================="
    echo "Meta: ${meta}"
    echo "Combined results: ${combined_results}"

    filter_doublets.py \\
        --input_h5ad ${raw_h5ad} \\
        --combined_results ${combined_results} \\
        --method ${params.doublet_consensus_method} \\
        --output_h5ad ${meta.id}_doublet_filtered.h5ad \\
        --summary_txt ${meta.id}_doublet_filter_summary.txt \\
        --image_prefix ${meta.id}_
    """
}
