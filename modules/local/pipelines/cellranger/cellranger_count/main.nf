process CR_PIPELINE {
    publishDir "${params.outdir}/CellRanger_pipeline/", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'


    container "quay.io/nf-core/cellranger:9.0.1"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path(cr_reference_dir)

    output:
    path("${meta.id}_count/outs"), emit: outs
    path "versions.yml",           emit: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba    TODO: add condition that singularity is not available
    // if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
    //     error "CELLRANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    // }

    // Index FASTQs, whether named by MERGE_FASTQS (<id>_merged_I1) or OAK_DEMUX (<id>_S1_L001_I1_001)
    def idx_files = [fastq_indices].flatten().findAll { it }
    def fastq_I1  = idx_files.find { it.name ==~ /.*_I1(_001)?\.f(ast)?q(\.gz)?/ }
    def fastq_I2  = idx_files.find { it.name ==~ /.*_I2(_001)?\.f(ast)?q(\.gz)?/ }
    def cr_prefix = "fastqs/${meta.id}_S1_L001"

    """
    echo "\n\n=============== CellRanger pipeline  ==============="
    echo "Sample ID: ${meta}"
    echo "Reference directory: ${cr_reference_dir}"

    # Cell Ranger only finds FASTQs named <sample>_S<n>_L<lane>_<read>_001.fastq.gz,
    # so link the inputs under those names, whatever they were called before
    mkdir -p fastqs
    ln -s "\$(pwd)/${fastq_BC_UMI}" ${cr_prefix}_R1_001.fastq.gz
    ln -s "\$(pwd)/${fastq_cDNA}" ${cr_prefix}_R2_001.fastq.gz
    ${fastq_I1 ? "ln -s \"\$(pwd)/${fastq_I1}\" ${cr_prefix}_I1_001.fastq.gz" : ''}
    ${fastq_I2 ? "ln -s \"\$(pwd)/${fastq_I2}\" ${cr_prefix}_I2_001.fastq.gz" : ''}
    ls -l fastqs/

    cellranger count \\
        --id=${meta.id}_count \\
        --transcriptome=${cr_reference_dir} \\
        --fastqs=fastqs \\
        --sample=${meta.id} \\
        --chemistry=auto \\
        --create-bam true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """
}
