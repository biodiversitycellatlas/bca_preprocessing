process MULTIQC {
    publishDir "${params.outdir}/summary_results", mode: 'copy'
    label 'process_single2'


    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/multiqc:1.30--d3e586af7b974fba"

    input:
    val(trigger)
    path(multiqc_config)

    output:
    path("*"),          emit: report
    path "versions.yml", emit: versions

    script:
    def config = multiqc_config ? "--config $multiqc_config" : ''

    """
    multiqc ${params.outdir} $config

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}
