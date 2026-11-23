process MERGE_REF_FASTA {
    publishDir path: { "${params.outdir}/genome" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_single2'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930' :
        'community.wave.seqera.io/library/python:3.14.2--0562a0df3245213a' }"

    input:
    path base_fasta
    path add_fasta

    output:
    path "ref.fasta",    emit: fasta
    path "versions.yml", emit: versions

    script:
    def do_merge = add_fasta ? true : false
    """
    if [ "$do_merge" = "true" ] && [ -s "${add_fasta}" ]; then
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
