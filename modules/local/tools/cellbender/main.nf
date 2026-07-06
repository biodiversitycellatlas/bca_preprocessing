process CELLBENDER {
    publishDir "${params.outdir}/cellbender/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label ( workflow.profile?.contains('gpu') ? 'process_gpu' : 'process_high' )
    debug true

    container "us.gcr.io/broad-dsde-methods/cellbender:latest"

    input:
    tuple val(meta), path(raw_h5ad)

    output:
    path("cellbender_output_report.html"),      emit: cb_html
    path("cellbender_output_metrics.csv"),      emit: cb_metrics
    path("cellbender_output.h5"),               emit: cb_output_h5
    path("cellbender_output_filtered.h5"),      emit: cb_output_filtered_h5
    path "versions.yml",                        emit: versions

    script:
    // Check if the workflow profile includes 'gpu' to determine whether to use GPU or CPU threads
    def use_gpu = task.ext.use_gpu ?: false
    def cuda_flag = use_gpu ? "--cuda" : "--cpu-threads ${task.cpus}"
    """
    echo "\n\n===============  Ambient RNA removal  ==============="
    echo "Sample ID: ${meta}"
    echo "Input files: ${raw_h5ad}"
    echo "Number of expected cells: ${meta.expected_cells}"
    echo "GPU Profile Active: ${use_gpu}"

    cellbender remove-background \\
        ${cuda_flag} \\
        --input ${raw_h5ad} \\
        --output cellbender_output.h5 \\
        --epochs 150 \\
        --expected-cells ${meta.expected_cells} \\
        --fpr 0.01 \\
        ${params.cellbender_extraargs ?: ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellbender: \$(python3 -c "import cellbender; print(cellbender.__version__)")
    END_VERSIONS
    """
}
