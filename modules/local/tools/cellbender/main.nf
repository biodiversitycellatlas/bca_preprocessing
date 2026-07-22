process CELLBENDER {
    publishDir "${params.outdir}/cellbender/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    debug true

    // Containers are the preferred way to run this module, and the only one that supports
    // GPU acceleration -- the bioconda build ships CPU-only torch. Container is selected
    // from task.ext.use_gpu, which conf/base.config sets on the process_gpu label: the
    // Broad image has CUDA-enabled torch, the Wave image is a slimmer CPU-only build.
    conda "${moduleDir}/environment.yml"
    container "${ task.ext.use_gpu ? 'us.gcr.io/broad-dsde-methods/cellbender:0.3.2' :
        workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/ebcf140f995f79fcad5c17783622e000550ff6f171771f9fc4233484ee6f63cf/data' :
        'community.wave.seqera.io/library/cellbender_webcolors:156d413fdfc16cdb' }"

    input:
    tuple val(meta), path(raw_h5ad)

    output:
    path("cellbender_output_report.html"),      emit: cb_html
    path("cellbender_output_metrics.csv"),      emit: cb_metrics
    path("cellbender_output.h5"),               emit: cb_output_h5
    path("cellbender_output_filtered.h5"),      emit: cb_output_filtered_h5
    path "versions.yml",                        emit: versions

    script:
    // GPU acceleration is enabled by running with `-profile gpu` (see conf/base.config).
    // --cpu-threads is always passed: it drives data loading, which stays on the CPU.
    def use_gpu = task.ext.use_gpu ?: false
    def cuda_flag = use_gpu ? "--cuda" : ""
    """
    echo "\n\n===============  Ambient RNA removal  ==============="
    echo "Sample ID: ${meta}"
    echo "Input files: ${raw_h5ad}"
    echo "Number of expected cells: ${meta.expected_cells}"
    echo "GPU acceleration: ${use_gpu ? 'enabled (--cuda)' : 'disabled (CPU only)'}"

    TMPDIR=. cellbender remove-background \\
        --cpu-threads ${task.cpus} \\
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
