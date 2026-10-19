process MERGE_REF_GTF {
    // The task's file stays ref.gtf, which STARSOLO_INDEX names its directory after; only the
    // published copy carries the suffix, so the GeneExt merge does not overwrite the standard one
    publishDir path: { "${params.outdir}/genome" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : (filename == 'ref.gtf' ? "ref${publish_suffix}.gtf" : filename) }
    label 'process_single2'

    input:
    path base_gtf
    path add_gtf
    val publish_suffix

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
