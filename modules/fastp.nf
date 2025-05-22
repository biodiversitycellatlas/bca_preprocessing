process FASTP {
    publishDir "${params.output_dir}/fastp", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)

    output:
    tuple val(meta), path("trimmed_${fastq_cDNA}"), path("trimmed_${fastq_BC_UMI}")

    script:
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
        --out2 trimmed_${fastq_BC_UMI}
    """
}
