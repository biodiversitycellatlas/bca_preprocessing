process CR_PIPELINE {
    publishDir "${params.resDir}/CellRanger_pipeline/\${basename_ref}", mode: 'copy', overwrite: false
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), path(fastq_files)
    path (cr_reference_dir)

    script:
    """
    echo "\n\n=============== CellRanger pipeline  ==============="
    echo "Sample ID: ${sample_id}"
    
    # Add CellRanger to PATH
    export PATH=${params.cellranger_dir}:$PATH

    # Define name of reference used during analysis
    basename_ref=\$(basename ${params.ref_star_gtf} .gtf)
    echo "basename: \${basename_ref}"

    cellranger count \\
        --id=${sample_id}_count \\
        --transcriptome=${cr_reference_dir} \\
        --fastqs=${params.resDir}/fastq/ \\
        --sample=${sample_id} \\
        --chemistry=auto \\
        --create-bam true
    """
}