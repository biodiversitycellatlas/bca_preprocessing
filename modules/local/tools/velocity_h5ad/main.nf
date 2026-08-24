process VELOCITY_H5AD {
    publishDir "${params.outdir}/anndata/${meta.id}/velocity", mode: 'copy'
    tag "${meta.id} | ${meta.mapping_method} | ${meta.datatype}"
    label 'process_single_mem2'

    // Memory tracks the size of the matrices being read, overrides process_single_mem2's flat
    // assignment. Coefficients live in params.dynamic_memory; remove the entry to fall back to
    // the plain label.
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.VELOCITY_H5AD, [matrix_dir], task.attempt, 12) }

    container 'oras://community.wave.seqera.io/library/scanpy:1.12--45f1dccaf83880df'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(matrix_dir)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml",             emit: versions

    script:
    // The two mappers store the three blocks differently, so the script is told which layout to read
    def source_arg = meta.mapping_method == 'alevin-fry' ? '--alevin-dir' : '--starsolo-dir'
    """
    velocity_matrices_to_h5ad.py \\
        ${source_arg} ${matrix_dir} \\
        --sample-id ${meta.id} \\
        --out ${meta.id}_${meta.datatype}_velocity.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        anndata: \$(python3 -c "import anndata; print(anndata.__version__)")
    END_VERSIONS
    """
}
