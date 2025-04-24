process FASTQC {
    publishDir "${params.output_dir}/fastqc", mode: 'copy'
    tag "${fastq_files}"

    input:
    tuple val(meta), path(fastq_files)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    echo "\n\n==================  FASTQC  =================="
    echo "Running FASTQC for ${fastq_files}"
    echo "Path: ${fastq_files}"

    fastqc ${fastq_files} 
    """
}
