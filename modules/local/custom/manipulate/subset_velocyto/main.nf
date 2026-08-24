process SUBSET_VELOCYTO_MATRICES {
    publishDir "${params.outdir}/mapping_STARsolo/${meta.id}/${meta.id}_Solo.out/Velocyto/filtered_secondderiv", mode: 'copy'
    label 'process_single_mem2'
    tag "${meta.id}"

    // Memory tracks the size of the three matrices being read, overrides process_single_mem2's
    // flat assignment. Coefficients live in params.dynamic_memory; remove the entry to fall
    // back to the plain label.
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.SUBSET_VELOCYTO_MATRICES, [velocyto_raw_dir], task.attempt, 10) }

    conda "${moduleDir}/environment.yml"

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
