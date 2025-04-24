process CR_PIPELINE {
    publishDir "${params.output_dir}/CellRanger_pipeline/${meta.id}", mode: 'copy', overwrite: false
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(fastqs)
    path(cr_reference_dir)

    script:
    """
    echo "\n\n=============== CellRanger pipeline  ==============="
    echo "Sample ID: ${meta}"
    echo "Reference directory: ${cr_reference_dir}"
    
    # Add CellRanger to PATH
    export PATH=${params.cellranger_dir}:$PATH

    clean_sample_id=\$(echo ${meta.id} | sed 's/_S1//')
    echo "cleaned sample id: \${clean_sample_id}"

    cellranger count \\
        --id=\${clean_sample_id}_count \\
        --transcriptome=${cr_reference_dir} \\
        --fastqs=${fastqs} \\
        --sample=\${clean_sample_id} \\
        --chemistry=auto \\
        --create-bam true
    """
}