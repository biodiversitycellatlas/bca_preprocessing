// =============  REF GENOME PARSE PIPELINE  ============= \\ 
// Creates the Genome Indeces required for mapping.        \\
// parse_refgenome: contains the reference genome for the  \\
// Parse Biosciences (split-pipe) pipeline.                \\

process REFGEN_PARSEBIO {
    output:
    path("*")

    script:
    """
    echo "\n\n==================  REF GENOME PARSE PIPELINE  =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    split-pipe -m mkref \\
        --genome_name ${params.species} \\
        --genes ${params.ref_parse_gtf} \\
        --fasta ${params.ref_fasta} \\
        --output_dir .
    """
}
