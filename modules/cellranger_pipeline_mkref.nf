process CR_PIPELINE_MKREF {
    publishDir "${params.output_dir}/genome/cellranger_ref", mode: 'copy'
    debug true

    output:
    path "*"

    script:
    """
    echo "\n\n=============== CellRanger pipeline - create REF  ==============="
    
    # Add CellRanger to PATH
    export PATH=${params.cellranger_dir}:$PATH


    cellranger mkref \\
        --fasta=${params.ref_fasta} \\
        --genes=${params.ref_gtf}
    """
}