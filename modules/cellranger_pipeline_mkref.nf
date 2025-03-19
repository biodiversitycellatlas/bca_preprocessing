process CR_PIPELINE_MKREF {
    publishDir "${params.resDir}/genome/cellranger_ref/\${basename_ref}", mode: 'copy'
    debug true

    output:
    path "${params.resDir}/genome/cellranger_ref/\${basename_ref}"

    script:
    """
    echo "\n\n=============== CellRanger pipeline - create REF  ==============="
    
    # Add CellRanger to PATH
    export PATH=${params.cellranger_dir}:$PATH

    # Define name of reference genome
    basename_ref=\$(basename ${params.ref_star_gtf})

    cellranger mkref \\
        --genome=\${basename_ref} \\
        --fasta=${params.ref_fasta} \\
        --genes=${params.ref_star_gtf}
    """
}