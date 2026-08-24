process COLLAPSE_ALEVIN_USA {
    publishDir "${params.outdir}/mapping_alevin/${meta.id}/gene_level_matrix/${meta.datatype}", mode: 'copy'
    tag "${meta.id} | ${meta.datatype} | ${usa_counts}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(alevin_mtx_dir)
    val usa_counts

    output:
    tuple val(meta), path("gene_level/"), emit: matrix
    path "versions.yml",                  emit: versions

    script:
    """
    collapse_alevin_usa.py \\
        --dir ${alevin_mtx_dir} \\
        --outdir gene_level/ \\
        --counts ${usa_counts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
