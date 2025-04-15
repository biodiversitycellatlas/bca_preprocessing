process CR_PIPELINE {
    publishDir "${params.output_dir}/CellRanger_pipeline/${sample_id}", mode: 'copy', overwrite: false
    tag "${sample_id}"
    debug true

    input:
    val(sample_id)
    path(cr_reference_dir)

    script:
    """
    echo "\n\n=============== CellRanger pipeline  ==============="
    echo "Sample ID: ${sample_id}"
    
    # Add CellRanger to PATH
    export PATH=${params.cellranger_dir}:$PATH

    clean_sample_id=\$(echo "$sample_id" | sed 's/_S1//')
    echo "cleaned sample id: \${clean_sample_id}"

    cellranger count \\
        --id=\${clean_sample_id}_count \\
        --transcriptome=${cr_reference_dir} \\
        --fastqs=${params.fastq_dir} \\
        --sample=\${clean_sample_id} \\
        --chemistry=auto \\
        --create-bam true
    """
}