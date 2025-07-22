process CR_PIPELINE_MKREF {
    publishDir "${params.output_dir}/genome/cellranger_ref", mode: 'copy'
    debug true

    output:
    path "*"

    script:
    """
    echo "\n\n=============== CellRanger pipeline - create REF  ==============="
    
    # Add CellRanger to PATH
    export PATH=${params.external_pipeline}:$PATH

    GTF_FILE="${params.ref_gtf_alt ?: params.ref_gtf}"

    cellranger mkref \\
        --genome cellranger_ref \\
        --fasta=${params.ref_fasta} \\
        --genes=\$GTF_FILE
    """
}