process SCALEBIO_PREPARE_BARCODES {
    publishDir path: { "${params.outdir}/demultiplex/scalebio/barcodes" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930' :
        'community.wave.seqera.io/library/python:3.14.2--0562a0df3245213a' }"

    input:
    // The user's barcode lists, one per segment in the layout's segment order
    path(barcode_lists, stageAs: 'lists/?/*')
    val(layout)

    output:
    path("barcodes.json"), emit: barcodes
    path("whitelist_*.txt"), emit: whitelists
    path "versions.yml", emit: versions

    script:
    def lists = (barcode_lists instanceof List ? barcode_lists : [barcode_lists]).join(' ')
    """
    echo -e "\\n\\n==================  ScaleBio: barcode lists  =================="
    scalebio_demux.py prepare \\
        --layout ${layout} \\
        --whitelists ${lists} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
