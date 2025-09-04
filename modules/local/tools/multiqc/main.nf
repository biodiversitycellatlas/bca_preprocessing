process MULTIQC {
    publishDir "${params.outdir}/summary_results", mode: 'copy'
    label 'process_single'
    

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.30--pyhdfd78af_1' :
        'biocontainers/multiqc:1.30--pyhdfd78af_1' }"

    input:
    val(trigger)
    path(multiqc_config)

    output:
    path("*")

    script:
    def config = multiqc_config ? "--config $multiqc_config" : ''

    """
    multiqc ${params.outdir} $config
    """
}
