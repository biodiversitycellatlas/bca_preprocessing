process MERGE_REF_FASTA {
    publishDir "${params.outdir}/genome", mode: 'copy'
    label 'process_single'

    input:
    path base_fasta
    path add_fasta

    output:
    path "ref.fasta",    emit: fasta
    path "versions.yml", emit: versions

    script:
    def do_merge = add_fasta ? true : false
    """
    if [ "$do_merge" = "true" ]; then
        cat ${base_fasta} ${add_fasta} > ref.fasta
    else
        cp ${base_fasta} ref.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/^GNU bash, version //; s/ .*\$//')
    END_VERSIONS
    """
}
