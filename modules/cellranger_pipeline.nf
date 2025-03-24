process CR_PIPELINE {
    publishDir "${params.resDir}/CellRanger_pipeline/\${basename_ref}", mode: 'copy', overwrite: false
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

    # Define name of reference used during analysis
    basename_ref=\$(basename ${params.ref_star_gtf} .gtf)
    echo "basename: \${basename_ref}"

    clean_sample_id=\$(echo "$sample_id" | sed 's/_S1//')
    echo "cleaned sample id: \${clean_sample_id}"

    cellranger count \\
        --id=\${clean_sample_id}_count \\
        --transcriptome=${cr_reference_dir} \\
        --fastqs=${params.resDir}/fastq/ \\
        --sample=\${clean_sample_id} \\
        --chemistry=auto \\
        --create-bam true
    """
}