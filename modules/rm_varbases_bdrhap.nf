process RM_VARBASES {
    tag "${meta.id}"
    debug true
    
    input:
    tuple val(meta), path(fastq_files)

    output:
    tuple val(meta), path("noVB_${meta.id}_R{1,2}_001.fastq.gz"), emit: fastq_noVB_files

    script:
    """
    # Activate the virtual environment
    source ~/cutadapt/bin/activate

    # Run cutadapt to remove the variable bases
    cutadapt \\
        -g "^A" -g "^GT" -g "^TCA" \\
        --quality-cutoff 0 \\
        -o noVB_${meta.id}_R1_001.fastq.gz -p noVB_${meta.id}_R2_001.fastq.gz \\
        ${meta.id}_R1_001.fastq.gz ${meta.id}_R2_001.fastq.gz

    # Deactivate the virtual environment after completion
    deactivate
    """
}