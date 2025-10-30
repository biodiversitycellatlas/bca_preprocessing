process SALMON_INDEX {
    publishDir "${params.outdir}/genome/salmon_index", mode: 'copy'
    label 'process_high'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.19.4--ha6fb395_0':
        'biocontainers/simpleaf:0.19.4--ha6fb395_0' }"

    input:
    path(splici_index_reference)

    output:
    path("salmon_index/")

    script:
    """
    echo "\n\n==================  SALMON INDEX =================="
    # Build reference index
    salmon index \\
        -t ${splici_index_reference}/*.fa \\
        -i ./salmon_index \\
        -k 31
    """
}
