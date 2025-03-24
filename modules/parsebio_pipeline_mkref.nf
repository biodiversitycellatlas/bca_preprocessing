process PARSEBIO_PIPELINE_MKREF {
    output:
    path("*")

    script:
    """
    echo "\n\n==================  REF GENOME PARSE PIPELINE  =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    split-pipe -m mkref \\
        --genome_name ref_splitpipe \\
        --genes ${params.ref_parse_gtf} \\
        --fasta ${params.ref_fasta} \\
        --output_dir .
    """
}
