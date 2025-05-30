process PARSEBIO_PIPELINE_MKREF {
    output:
    path("*")

    script:
    """
    echo "\n\n==================  REF GENOME PARSE PIPELINE  =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    # Check if a specific reference GTF file is provided, otherwise use the default
    GTF_FILE="${params.ref_gtf_alt ?: params.ref_gtf}"

    # Create reference genome index files
    split-pipe -m mkref \\
        --genome_name ref_splitpipe \\
        --genes \$GTF_FILE \\
        --fasta ${params.ref_fasta} \\
        --output_dir .
    """
}
