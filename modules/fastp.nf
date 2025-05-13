process FASTP {
    publishDir "${params.output_dir}/fastp/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.trim.fastq.gz")

    script:
    """
    echo "\n\n==================  TRIM FASTQs WITH FASTP  =================="
    echo "Sample ID: ${meta}"
    echo "Processing files: ${fastqs}"

    fastp \\
        --html fastp_${meta.id}.html \\
        --json fastp_${meta.id}.json \\
        --thread 8 \\
        --in1 ${fastqs[0]} \\
        --in2 ${fastqs[1]} \\
        --out1 ${meta.id}_R1.trim.fastq.gz \\
        --out2 ${meta.id}_R2.trim.fastq.gz
    """
}
