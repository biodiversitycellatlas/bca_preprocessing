process RM_VARBASES {
    tag "${meta.id}"
    label 'process_medium'


    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/cutadapt:5.2--5505472a18e9cce0"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    tuple val(meta), path("noVB_${meta.id}_R2_001.fastq.gz"), path("noVB_${meta.id}_R1_001.fastq.gz"), path(fastq_indices), path(input_file)

    script:
    """
    # Run cutadapt to remove the variable bases
    cutadapt \\
        -g "^A" -g "^GT" -g "^TCA" \\
        --quality-cutoff 0 \\
        -o noVB_${meta.id}_R1_001.fastq.gz -p noVB_${meta.id}_R2_001.fastq.gz \\
        ${fastq_BC_UMI} ${fastq_cDNA}
    """
}
