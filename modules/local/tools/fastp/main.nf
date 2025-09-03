process FASTP {
    publishDir "${params.outdir}/fastp", mode: 'copy'
    tag "${meta.id}"
    label 'process_low'
    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/88/889a182b8066804f4799f3808a5813ad601381a8a0e3baa4ab8d73e739b97001/data' :
        'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690' }"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)

    output:
    tuple val(meta), path("trimmed_${fastq_cDNA}"), path("trimmed_${fastq_BC_UMI}")

    script:
    // Retrieve fastp settings from conf/seqtech_parameters.config
    def fastp_settings = params.seqtech_parameters[params.protocol].fastp
    def fastp_args = fastp_settings instanceof List ? fastp_settings.join(' ') : fastp_settings
    
    """
    echo "\n\n==================  TRIM FASTQs WITH FASTP  =================="
    echo "Metadata: ${meta}"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"

    fastp \\
        --html fastp_${meta.id}.html \\
        --json fastp_${meta.id}.json \\
        --thread 8 \\
        --in1 ${fastq_cDNA} \\
        --in2 ${fastq_BC_UMI} \\
        --out1 trimmed_${fastq_cDNA} \\
        --out2 trimmed_${fastq_BC_UMI} \\
        ${fastp_args}
    """
}
