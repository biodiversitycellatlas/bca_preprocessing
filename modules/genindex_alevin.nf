process GENINDEX_ALEVIN { 
    publishDir "${params.resDir}/genome/genindex_alevin", mode: 'copy'
    debug true

    output:
    path("*")
       
    script:
    """
    echo "\n\n==================  GENOME INDEX ALEVIN =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Reference fasta: ${params.ref_fasta}"
    echo "Reference ref_star_gtf: ${params.ref_star_gtf}"

    # Retrieve the first accession number
    first_fastq=\$(ls "${params.fastq_dir}/" | head -n1)  

    # Calculate read length using the first read from the first fastq file
    readlen=\$(zcat ${params.fastq_dir}/\${first_fastq} | awk 'NR==2 {print length(\$0)}') 

    echo "First FASTQ file: \${first_fastq}"
    echo "Read length: \${readlen}"

    # Create splici reference
    Rscript ${params.baseDir}/bin/salmon_create_splici_ref.R \\
        --ref_fasta ${params.ref_fasta} \\
        --ref_gtf ${params.ref_star_gtf} \\
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
