process GENINDEX_ALEVIN {   
    input:
    path(ref_star_gtf)

    output:
    path("index")                   , emit: index
    path("ref/*_t2g_3col.tsv")      , emit: transcript_tsv
    path("*")                       , emit: salmon
       
    script:
    """
    echo "\n\n==================  GENOME INDEX ALEVIN =================="
    echo "Reference fasta: ${params.ref_fasta}"
    
    # simpleaf configuration
    simpleaf set-paths

    # Build reference index
    simpleaf index \\
        --fasta ${params.ref_fasta} \\
        --gtf \${ref_star_gtf}
    """
}
