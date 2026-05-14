process MERGE_REF_GTF {
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path base_gtf
    path add_gtf

    output:
    path "ref.gtf"

    script:
    def do_merge = add_gtf ? true : false
    """
    if [ "$do_merge" = "true" ]; then
        cat ${base_gtf} ${add_gtf} > ref.gtf
    else
        cp ${base_gtf} ref.gtf
    fi
    """
}
