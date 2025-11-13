process RM_VARBASES {
    tag "${meta.id}"
    label 'process_medium'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/17/1758869538eb8e658077cc14cd7a4e76fd9b6d73d3a68f85a70bf292e39e27c5/data' :
        'community.wave.seqera.io/library/cutadapt:5.0--991bbd2e184b7014' }"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    tuple val(meta), path("noVB_${meta.id}_R2_001.fastq.gz"), path("noVB_${meta.id}_R1_001.fastq.gz"), path(input_file), emit: fastq_noVB_files

    script:
    """
    # Run cutadapt to remove the variable bases
    cutadapt \\
        -g "^A" -g "^GT" -g "^TCA" \\
        --quality-cutoff 0 \\
        -o noVB_${meta.id}_R1_001.fastq.gz -p noVB_${meta.id}_R2_001.fastq.gz \\
        ${meta.id}_R1_001.fastq.gz ${meta.id}_R2_001.fastq.gz
    """
}
