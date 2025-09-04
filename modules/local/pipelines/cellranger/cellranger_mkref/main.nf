process CR_PIPELINE_MKREF {
    publishDir "${params.outdir}/genome/cellranger_ref", mode: 'copy'
    label 'process_medium'
    

    container "quay.io/nf-core/cellranger:9.0.1"

    output:
    path "*"

    script:
    // Exit if running this module with -profile conda / -profile mamba  TODO: add condition that singularity is not available
    // if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
    //     error "CELLRANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    // }
    
    """
    echo "\n\n=============== CellRanger pipeline - create REF  ==============="

    GTF_FILE="${params.ref_gtf_alt ?: params.ref_gtf}"

    cellranger mkref \\
        --genome cellranger_ref \\
        --fasta=${params.ref_fasta} \\
        --genes=\$GTF_FILE
    """
}