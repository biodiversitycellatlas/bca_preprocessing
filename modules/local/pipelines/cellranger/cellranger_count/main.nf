process CR_PIPELINE {
    publishDir "${params.outdir}/CellRanger_pipeline/", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'
    

    container "quay.io/nf-core/cellranger:9.0.1"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(input_file)
    path(cr_reference_dir)

    output:
    path("${meta.id}_count/outs")

    script:
    // Exit if running this module with -profile conda / -profile mamba    TODO: add condition that singularity is not available
    // if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
    //     error "CELLRANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    // }

    """
    echo "\n\n=============== CellRanger pipeline  ==============="
    echo "Sample ID: ${meta}"
    echo "Reference directory: ${cr_reference_dir}"
    
    cellranger count \\
        --id=${meta.id}_count \\
        --transcriptome=${cr_reference_dir} \\
        --fastqs=. \\
        --sample=${meta.id} \\
        --chemistry=auto \\
        --create-bam true
    """
}