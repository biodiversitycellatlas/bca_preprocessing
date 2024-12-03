// =============  REF GENOME PARSE PIPELINE  ============= \\ 
// Creates the Genome Indeces required for mapping.        \\
// parse_refgenome: contains the reference genome for the  \\
// Parse Biosciences (split-pipe) pipeline.                \\

process REFGEN_PARSEBIO {
    publishDir "${params.resDir}/genome/parse_refgenome", mode: 'symlink'
    debug true

    output:
    path("*")

    script:
    """
    echo "\n\n==================  REF GENOME PARSE PIPELINE  =================="
    split-pipe -m mkref \\
        --genome_name ${params.species} \\
        --genes ${params.ref_parse_gtf} \\
        --fasta ${params.ref_fasta} \\
        --output_dir .
    """
}
