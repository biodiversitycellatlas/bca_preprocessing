process GENINDEX_ALEVIN { 
    publishDir "${params.output_dir}/genome/genindex_alevin", mode: 'copy'
    debug true

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)

    output:
    path("*")
       
    script:
    """
    echo "\n\n==================  GENOME INDEX ALEVIN =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Reference fasta: ${params.ref_fasta}"
    echo "Reference ref_gtf: ${params.ref_gtf}"

    # Calculate read length using the first read from the first fastq file
    readlen=\$(zcat ${fastq_cDNA} | awk 'NR==2 {print length(\$0)}') 
    echo "Read length: \${readlen}"

    # Create splici reference
    Rscript ${launchDir}/bin/salmon_create_splici_ref.R \\
        --ref_fasta ${params.ref_fasta} \\
        --ref_gtf ${params.ref_gtf} \\
        --readlen \${readlen} \\
        --flanklen 5 \\
        --prefix "transcriptome_splici" \\
        --out_dir ./splici_index_reference
    
    # Define the reference fasta file created by the R script
    ref_fasta=\$(ls ./splici_index_reference/*.fa)

    # Build reference index
    salmon index \\
        -t \${ref_fasta} \\
        -i ./salmon_index \\
        -p 32
    """
}
