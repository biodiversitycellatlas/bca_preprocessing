process SALMON_SPLICI {
    publishDir "${params.outdir}/genome/salmon_splici", mode: 'copy'
    label 'process_high'


    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    path("splici_index_reference/")

    script:
    """
    echo "\n\n==================  SALMON SPLICI =================="
    echo "Reference fasta: ${params.ref_fasta}"
    echo "Reference ref_gtf: ${params.ref_gtf}"

    # Calculate read length using the first read from the first fastq file
    readlen=\$(zcat ${fastq_cDNA} | awk 'NR==2 {print length(\$0)}')
    echo "Read length: \${readlen}"

    # Create splici reference
    salmon_create_splici_ref.R \\
        --ref_fasta ${params.ref_fasta} \\
        --ref_gtf ${params.ref_gtf} \\
        --readlen \${readlen} \\
        --flanklen 5 \\
        --prefix "transcriptome_splici" \\
        --out_dir ./splici_index_reference

    # Define the reference fasta file created by the R script
    # ref_fasta=\$(ls ./splici_index_reference/*.fa)
    """
}
