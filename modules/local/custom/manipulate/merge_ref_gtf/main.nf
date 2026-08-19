process MERGE_REF_GTF {
    publishDir "${params.outdir}/genome", mode: 'copy'
    label 'process_single'

    input:
    path base_gtf
    path add_gtf

    output:
    path "ref.gtf",      emit: gtf
    path "versions.yml", emit: versions

    script:
    def do_merge = add_gtf ? true : false
    """
    if [ "$do_merge" = "true" ]; then
        cat ${base_gtf} ${add_gtf} > ref.gtf
    else
        cp ${base_gtf} ref.gtf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/^GNU bash, version //; s/ .*\$//')
    END_VERSIONS
    """
}
