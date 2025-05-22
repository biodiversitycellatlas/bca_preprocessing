process FASTQC {
    publishDir "${params.output_dir}/fastqc", mode: 'copy'
    tag "${fastq_cDNA}, ${fastq_BC_UMI}"

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
