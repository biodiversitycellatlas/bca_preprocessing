process MTX_TO_10X {
    publishDir path: { "${params.outdir}/doublet_filtering/${meta.id}/10x_export" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${meta.mapping_method}"

    // Memory tracks the size of the matrix being read. Coefficients live in
    // params.dynamic_memory; without an entry the request is process_single_mem2's 10 GB.
    // Deliberately no resource label, since a withLabel setting would take precedence over
    // these: cpus and time are process_single_mem2's.
    cpus   1
    time   { 30.m * task.attempt }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.MTX_TO_10X, [mtx, barcodes, features], task.attempt, 10) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_pandas_scipy:9829b20b4522c6a0' :
        'community.wave.seqera.io/library/python_numpy_pandas_scipy:bb8d35a0908056f1' }"

    input:
    tuple val(meta), path(mtx), path(barcodes), path(features)

    output:
    tuple val(meta), path("${meta.id}_10x"), emit: tenx_dir

    script:
    """
    echo -e "\\n\\n==================  MTX to 10x export =================="
    echo "Meta: ${meta}"

    mtx_to_10x.py \\
        --mtx ${mtx} \\
        --barcodes ${barcodes} \\
        --features ${features} \\
        --outdir ${meta.id}_10x
    """
}
