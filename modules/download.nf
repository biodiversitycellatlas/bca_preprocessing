process DOWNLOAD_DATA {
    tag "${sample_id}"
    
    input:
    val sample_id

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_001.fastq.gz"), emit: fastq_files

    script:
    """
    echo "sample id: ${sample_id}"
    ln -s "${params.fastq_dir}/${sample_id}_R1_001.fastq.gz" .
    ln -s "${params.fastq_dir}/${sample_id}_R2_001.fastq.gz" .
    """
}
