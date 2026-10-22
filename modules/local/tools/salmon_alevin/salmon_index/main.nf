process SALMON_INDEX {
    publishDir path: { "${params.outdir}/genome" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_high'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/salmon:1.10.3--726401738a281398' :
        'community.wave.seqera.io/library/salmon:1.10.3--6ea6d56a4b0c26b8' }"

    input:
    path(splici_index_reference)
    val(index_suffix)

    output:
    path("salmon_index${index_suffix}/"), emit: salmon_index
    path "versions.yml",   emit: versions


    script:
    """
    echo -e "\\n\\n==================  SALMON INDEX =================="
    # Build reference index
    salmon index \\
        -t ${splici_index_reference}/*.fa \\
        -i ./salmon_index${index_suffix} \\
        -k 31

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}
