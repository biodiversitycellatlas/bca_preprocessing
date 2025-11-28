process SALMON_INDEX {
    publishDir "${params.outdir}/genome", mode: 'copy'
    label 'process_high'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.19.4--ha6fb395_0':
        'quay.io/biocontainers/simpleaf:0.19.4--ha6fb395_0' }"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path(splici_index_reference)

    output:
    path("salmon_index_${meta.id}/"), emit: salmon_index
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file), emit: salmon_samplesheet


    script:
    """
    echo "\n\n==================  SALMON INDEX =================="
    # Build reference index
    salmon index \\
        -t ${splici_index_reference}/*.fa \\
        -i ./salmon_index_${meta.id} \\
        -k 31
    """
}
