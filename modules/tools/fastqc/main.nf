process FASTQC {
    publishDir "${params.output_dir}/fastqc", mode: 'copy'
    tag "${fastq_cDNA}, ${fastq_BC_UMI}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc ${fastq_cDNA} 
    fastqc ${fastq_BC_UMI}
    """
}
