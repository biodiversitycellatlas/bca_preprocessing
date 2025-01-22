// ==================  FASTQC  ================== \\ 
// Generate Quality Control reports using FASTQC  \\

process FASTQC {
    publishDir "${params.resDir}/fastqc", mode: 'copy'
    tag "${fastq_files}"

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    path("*_fastqc.*")

    script:
    """
    echo "\n\n==================  FASTQC  =================="
    echo "Running FASTQC for ${fastq_files}"
    echo "Path: ${fastq_files}"

    fastqc ${fastq_files} 
    """
}
