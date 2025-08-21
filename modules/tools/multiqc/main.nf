process MULTIQC {
    publishDir "${params.output_dir}/summary_results", mode: 'copy'
    label 'process_single'
    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.30--pyhdfd78af_1' :
        'biocontainers/multiqc:1.30--pyhdfd78af_1' }"

    input:
    val(trigger)

    output:
    path("*")

    script:
    """
    multiqc ${params.output_dir} --config ${moduleDir}/bin/multiqc_config.yml
    """
}
