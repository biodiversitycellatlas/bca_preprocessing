process PARSEBIO_PIPELINE_MKREF {
    output:
    path("*")

    script:
    """
    echo "\n\n==================  REF GENOME PARSE PIPELINE  =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    # Install the required python packages if not already present
    pip install ${params.external_pipeline} --no-cache-dir

    GTF_FILE="${params.ref_gtf_alt ?: params.ref_gtf}"

    split-pipe -m mkref \\
        --genome_name ref_splitpipe \\
        --genes \$GTF_FILE \\
        --fasta ${params.ref_fasta} \\
        --output_dir .
    """
}
