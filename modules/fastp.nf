process FASTP {
    publishDir "${params.output_dir}/fastp/${sample_id}", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    path("*.trim.fastq.gz")

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    """
    echo "\n\n==================  TRIM FASTQs WITH FASTP  =================="
    echo "Sample ID: ${sample_id}"
    echo "Processing files: ${fastq_files}"

    fastp \\
        --html fastp_${sample_id}.html \\
        --json fastp_${sample_id}.json \\
        --thread 8 \\
        --in1 ${r1_fastq} \\
        --in2 ${r2_fastq} \\
        --out1 ${sample_id}_R1.trim.fastq.gz \\
        --out2 ${sample_id}_R2.trim.fastq.gz \\
        >& {log}
    """
}
