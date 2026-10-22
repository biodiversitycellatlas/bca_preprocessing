process DOUBLET_FILTER {
    publishDir path: { "${params.outdir}/doublet_filtering/${meta.id}/${meta.datatype}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${meta.mapping_method} | ${meta.datatype}"

    // Memory tracks the size of the matrix being sliced. Coefficients live in
    // params.dynamic_memory; without an entry the request is process_low's 12 GB.
    // Deliberately no resource label, since a withLabel setting would take precedence over
    // these: cpus and time are process_low's.
    cpus   { 2 * task.attempt }
    time   { 4.h * task.attempt }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.DOUBLET_FILTER, [mtx, barcodes, features], task.attempt, 12) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_scipy_matplotlib:3ab628a840bacec5' :
        'community.wave.seqera.io/library/python_pandas_numpy_scipy_matplotlib:ef1cecd8d63c13b2' }"

    input:
    tuple val(meta), path(mtx), path(barcodes), path(features), path(combined_results)

    output:
    tuple val(meta), path("${meta.id}_${meta.datatype}_doublet_filtered"),      emit: matrix
    path("${meta.id}_${meta.datatype}_doublet_filter_summary.txt"),             emit: summary
    path("${meta.id}_${meta.datatype}_doublet_filtering_summary.png"),          emit: filtering_summary_plot

    script:
    """
    echo -e "\\n\\n==================  Doublet consensus filter =================="
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
