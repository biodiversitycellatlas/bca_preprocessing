process APPLY_METACELL_FILTER {
    publishDir path: { "${params.outdir}/metacell_filtering/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${meta.mapping_method}"
    label 'process_single_mem2'

    container 'oras://community.wave.seqera.io/library/scanpy:1.12--45f1dccaf83880df'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cells_h5ad)
    path(selection)

    output:
    tuple val(meta), path("${meta.id}_final"),                                    emit: matrix_dir
    tuple val(meta), path("${meta.id}_final/${meta.id}_final.h5ad"),              emit: h5ad
    tuple val(meta), path("${meta.id}_final/${meta.id}_filter_summary.json"),     emit: summary
    path "versions.yml",                                                          emit: versions

    script:
    """
    apply_metacell_filter.py \\
        --cells_h5ad ${cells_h5ad} \\
        --selection ${selection} \\
        --sample_id ${meta.id} \\
        --outdir ${meta.id}_final \\
        --gene_mode ${params.metacell_gene_mode ?: 'drop'} \\
        --compression ${params.h5ad_compression ?: 'gzip'} \\
        --versions_yml versions.yml \\
        --process_name "${task.process}"
    """

    stub:
    """
    mkdir -p ${meta.id}_final
    touch ${meta.id}_final/matrix.mtx ${meta.id}_final/barcodes.tsv ${meta.id}_final/features.tsv
    touch ${meta.id}_final/${meta.id}_final.h5ad ${meta.id}_final/${meta.id}_final_metacells.h5ad
    echo '{}' > ${meta.id}_final/${meta.id}_filter_summary.json
    touch versions.yml
    """
}
