process CR_PIPELINE {
    publishDir "${params.output_dir}/CellRanger_pipeline/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)
    path(cr_reference_dir)

    script:
    """
    echo "\n\n=============== CellRanger pipeline  ==============="
    echo "Sample ID: ${meta}"
    echo "Reference directory: ${cr_reference_dir}"
    
    # Add CellRanger to PATH
    export PATH=${params.cellranger_dir}:$PATH

    # Derive the parent directory of the first FASTQ
    fastq_parent=\$(dirname ${fastq_cDNA})
    echo "Using FASTQ directory: \$fastq_parent"
    
    cellranger count \\
        --id=${meta.id}_count \\
        --transcriptome=${cr_reference_dir} \\
        --fastqs=\$fastq_parent \\
        --sample=${meta.id} \\
        --chemistry=auto \\
        --create-bam true
    """
}