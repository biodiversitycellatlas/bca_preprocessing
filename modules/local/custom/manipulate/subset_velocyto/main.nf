process SUBSET_VELOCYTO_MATRICES {
    publishDir path: { "${params.outdir}/mapping_STARsolo/${meta.id}/${meta.id}_Solo.out/Velocyto/filtered_secondderiv" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"

    // Memory tracks the size of the three matrices being read. Coefficients live in
    // params.dynamic_memory; without an entry the request is process_single_mem2's 10 GB.
    // Deliberately no resource label, since a withLabel setting would take precedence over
    // these: cpus and time are process_single_mem2's.
    cpus   1
    time   { 30.m * task.attempt }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.SUBSET_VELOCYTO_MATRICES, [velocyto_raw_dir], task.attempt, 10) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_pandas_scipy:eeb34b63825ecb6a' :
        'community.wave.seqera.io/library/python_numpy_pandas_scipy:b12b357472dfbd06' }"

    input:
    tuple val(meta), path(velocyto_raw_dir), path(filtered_barcodes)

    output:
    tuple val(meta), path("filtered/"), emit: filtered_matrix
    path "versions.yml",                emit: versions

    script:
    """
    echo "Subsetting the Velocyto matrices of ${meta.id} to \$(wc -l < ${filtered_barcodes}) called cells"

    subset_matrices_to_cells.py \\
        --dir ${velocyto_raw_dir} \\
        --barcodes ${filtered_barcodes} \\
        --outdir filtered/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
