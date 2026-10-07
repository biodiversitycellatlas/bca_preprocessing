process MULTIQC {
    publishDir path: { "${params.outdir}/summary_results" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_single2'


    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/multiqc:1.30--d3e586af7b974fba"

    input:
    path(multiqc_files, stageAs: "?/*")
    tuple val(salmon_ids), path(salmon_meta, stageAs: "salmon_meta/?/meta_info.json")
    path(multiqc_config)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"),        emit: data
    path("multiqc_plots"),       emit: plots, optional: true
    path "versions.yml",         emit: versions

    script:
    def config = multiqc_config ? "--config $multiqc_config" : ''

    // Rebuild '<id>_run/aux_info/meta_info.json', where MultiQC reads the sample name from
    def salmon_list    = salmon_meta instanceof List ? salmon_meta : [salmon_meta]
    def salmon_restage = [salmon_ids, salmon_list].transpose().collect { id, meta_info ->
        "mkdir -p salmon/${id}_run/aux_info && cp -L ${meta_info} salmon/${id}_run/aux_info/meta_info.json"
    }.join('\n    ')

    """
    ${salmon_restage}

    multiqc --force $config --ignore salmon_meta .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}
