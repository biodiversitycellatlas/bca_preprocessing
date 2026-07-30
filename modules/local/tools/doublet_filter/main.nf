process DOUBLET_FILTER {
    publishDir "${params.outdir}/doublet_filtering/${meta.id}/${meta.datatype}", mode: 'copy'
    tag "${meta.id} | ${meta.mapping_method} | ${meta.datatype}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input_h5ad), path(combined_results)

    output:
    tuple val(meta), path("${meta.id}_${meta.datatype}_doublet_filtered.h5ad"), emit: h5ad
    path("${meta.id}_${meta.datatype}_doublet_filter_summary.txt"),             emit: summary
    path("${meta.id}_${meta.datatype}_doublet_filtering_summary.png"),          emit: filtering_summary_plot

    script:
    """
    echo "\n\n==================  Doublet consensus filter =================="
    echo "Meta: ${meta}"
    echo "Input h5ad: ${input_h5ad}"
    echo "Combined results: ${combined_results}"

    filter_doublets.py \\
        --input_h5ad ${input_h5ad} \\
        --combined_results ${combined_results} \\
        --method ${params.doublet_consensus_method} \\
        --output_h5ad ${meta.id}_${meta.datatype}_doublet_filtered.h5ad \\
        --summary_txt ${meta.id}_${meta.datatype}_doublet_filter_summary.txt \\
        --image_prefix ${meta.id}_${meta.datatype}_
    """
}
