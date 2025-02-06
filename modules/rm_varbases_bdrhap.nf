// =============  REMOVE VARIABLE BASES FROM FASTQ ============= \\ 

process RM_VARBASES {
    tag "${sample_id}"
    debug true
    
    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val(sample_id), path("noVB_${sample_id}_R{1,2}_001.fastq.gz"), emit: fastq_noVB_files

    script:
    """
    # Activate the virtual environment
    source ~/cutadapt/bin/activate

    # Run cutadapt to remove the variable bases
    cutadapt \\
        -g "^A" -g "^GT" -g "^TCA" \\
        --quality-cutoff 0 \\
        -o noVB_${sample_id}_R1_001.fastq.gz -p noVB_${sample_id}_R2_001.fastq.gz \\
        ${sample_id}_R1_001.fastq.gz ${sample_id}_R2_001.fastq.gz

    # Deactivate the virtual environment after completion
    deactivate
    """
}