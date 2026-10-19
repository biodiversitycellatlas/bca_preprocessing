process MC_GENE_TABLE {
    publishDir path: { "${params.outdir}/metacells/gene_table" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_single'

    container "oras://community.wave.seqera.io/library/numpy_pandas_python:b78d8279027c3c38"
    conda "${moduleDir}/environment.yml"

    input:
    path(ref_gtf)
    path(gene_annotation)

    output:
    path("gene_table.tsv"),        emit: gene_table
    path("gene_table_stats.json"), emit: stats
    path "versions.yml",           emit: versions

    script:
    def annotation_arg = gene_annotation
        ? "--annotation ${gene_annotation} --annotation_id_col '${params.gene_annotation_id_col}' --annotation_pfam_col '${params.gene_annotation_pfam_col}'"
        : ""
    """
    build_gene_table.py \\
        --gtf ${ref_gtf} \\
        --mt_contig ${params.mt_contig ?: ''} \\
        --rrna_pattern '${params.grep_rrna ?: 'rRNA'}' \\
        ${annotation_arg} \\
        --output gene_table.tsv \\
        --stats gene_table_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    printf 'gene_id\\tgene_name\\tchrom\\tbiotype\\tis_mito\\tis_rrna\\tpfam\\n' > gene_table.tsv
    echo '{}' > gene_table_stats.json
    touch versions.yml
    """
}
