process DOWNLOAD_DATA {
    tag { meta.id }
    
    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("${meta.id}_S1_R{1,2}_001.fastq.gz")

    script:
    """
    echo "Metadata: ${meta}"
    echo "Fastq files: ${fastqs}"
    
    ln -s ${fastqs[0]} ${meta.id}_S1_R1_001.fastq.gz
    ln -s ${fastqs[1]} ${meta.id}_S1_R2_001.fastq.gz
    """
}
