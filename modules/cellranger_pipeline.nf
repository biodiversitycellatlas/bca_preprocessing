process CR_PIPELINE {
    publishDir "${params.output_dir}/CellRanger_pipeline/", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)
    path(cr_reference_dir)

    output:
    path("${meta.id}_count")

    script:
    """
    echo "\n\n=============== CellRanger pipeline  ==============="
    echo "Sample ID: ${meta}"
    echo "Reference directory: ${cr_reference_dir}"
    
    # Add CellRanger to PATH
    export PATH=${params.cellranger_dir}:$PATH

    # Create a temporary directory for symlinked FASTQ files
    mkdir -p fastq_temp

    # Extract base names (without .fastq.gz extension)
    base_bc_umi=\$(basename ${fastq_BC_UMI} .fastq.gz)
    base_cdna=\$(basename ${fastq_cDNA} .fastq.gz)

    # Create Cell Ranger-compliant symlinks
    ln -s \$(realpath ${fastq_BC_UMI}) fastq_temp/\${base_bc_umi}_S0_R1_001.fastq.gz
    ln -s \$(realpath ${fastq_cDNA})  fastq_temp/\${base_cdna}_S0_R2_001.fastq.gz

    echo "Symbolic links created in fastq_temp:"
    ls -lh fastq_temp
    
    cellranger count \\
        --id=${meta.id}_count \\
        --transcriptome=${cr_reference_dir} \\
        --fastqs=fastq_temp \\
        --sample=${meta.id} \\
        --chemistry=auto \\
        --create-bam true
    """
}