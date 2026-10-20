process MTX_TO_H5AD {
    publishDir path: { "${params.outdir}/anndata/${meta.id}/${meta.datatype}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${meta.mapping_method} | ${meta.datatype}"

    // Memory tracks the size of the matrix being read -- both of them, where CellSweep's
    // denoised counts are merged in. Coefficients live in params.dynamic_memory; without an
    // entry the request is process_single_mem2's 10 GB. Deliberately no resource label,
    // since a withLabel setting would take precedence over these: cpus and time are
    // process_single_mem2's.
    cpus   1
    time   { 30.m * task.attempt }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.MTX_TO_H5AD, [mtx, barcodes, features, cellsweep_h5ad], task.attempt, 10) }

    container 'oras://community.wave.seqera.io/library/scanpy:1.12--45f1dccaf83880df'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(mtx), path(barcodes), path(features), path(doublet_results), path(cellsweep_h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml",             emit: versions

    script:
    def doublet_arg   = doublet_results ? "--doublet_results ${doublet_results} --doublet_method ${params.doublet_consensus_method}" : ""
    def cellsweep_arg = cellsweep_h5ad  ? "--cellsweep_h5ad ${cellsweep_h5ad}" : ""
    """
    echo -e "\\n\\n==================  MTX to h5ad =================="
    echo "Meta: ${meta}"
    echo "Matrix: ${mtx}"
    echo "Doublet results: ${doublet_results ?: 'none'}"
    echo "CellSweep results: ${cellsweep_h5ad ?: 'none'}"

    # scanpy's numba functions cache next to their source, which a container image keeps read-only
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache

    mtx_to_h5ad.py \\
        --mtx ${mtx} \\
        --barcodes ${barcodes} \\
        --features ${features} \\
        --sample_id ${meta.id} \\
        --output_h5ad ${meta.id}_${meta.datatype}.h5ad \\
        --compression ${params.h5ad_compression ?: 'gzip'} \\
        --versions_yml versions.yml \\
        --process_name "${task.process}" \\
        ${doublet_arg} \\
        ${cellsweep_arg}
    """
}
