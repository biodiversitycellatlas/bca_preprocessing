process CELLBENDER {
    publishDir "${params.outdir}/cellbender/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label ( workflow.profile?.contains('gpu') ? 'process_gpu' : 'process_high' )
    debug true

    container "us.gcr.io/broad-dsde-methods/cellbender:latest"

    input:
    tuple val(meta), path(starsolo_genefull50_raw)

    output:
    path("cellbender_output*")

    script:
    // Check if the workflow profile includes 'gpu' to determine whether to use GPU or CPU threads
    def use_gpu = task.ext.use_gpu ?: false
    def cuda_flag = use_gpu ? "--cuda" : "--cpu-threads ${task.cpus}"
    """
    echo "\n\n===============  Ambient RNA removal  ==============="
    echo "Sample ID: ${meta}"
    echo "STARsolo files: ${starsolo_genefull50_raw}"
    echo "Number of expected cells: ${meta.expected_cells}"
    echo "GPU Profile Active: ${use_gpu}"

    # Copy features file as cellbender expects the file to be named genes.tsv
    cp ${starsolo_genefull50_raw}/features.tsv ${starsolo_genefull50_raw}/genes.tsv

    cellbender remove-background \\
        ${cuda_flag} \\
        --input ${starsolo_genefull50_raw} \\
        --output cellbender_output.h5 \\
        --epochs 150 \\
        --expected-cells ${meta.expected_cells} \\
        --fpr 0.01 \\
        ${params.cellbender_extraargs ?: ''}
    """
}
