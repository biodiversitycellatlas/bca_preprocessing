process GENINDEX_ALEVIN {   
    output:
    path("index")                   , emit: index
    path("ref/*_t2g_3col.tsv")      , emit: transcript_tsv
    path("*")                       , emit: salmon
       
    script:
    """
    echo "\n\n==================  GENOME INDEX ALEVIN =================="
    echo "Reference fasta: ${params.ref_fasta}"
    echo "Reference ref_star_gtf: ${params.ref_star_gtf}"
    
    # Filter transcripts spanning multiple chromosomes or strands
    awk '\$3 == "exon" {print \$12, \$1}' ${params.ref_star_gtf} | sort | uniq | awk '{print \$1}' | uniq -c | awk '\$1 > 1 {print \$2}' > transcripts_multiple_chromosomes.txt
    awk '\$3 == "exon" {print \$12, \$7}' ${params.ref_star_gtf} | sort | uniq | awk '{print \$1}' | uniq -c | awk '\$1 > 1 {print \$2}' > transcripts_multiple_strands.txt
    awk '\$3 == "exon" {print \$10, \$1}' ${params.ref_star_gtf} | sort | uniq | awk '{print \$1}' | uniq -c | awk '\$1 > 1 {print \$2}' > genes_multiple_chromosomes.txt
    awk '\$3 == "exon" {print \$10, \$7}' ${params.ref_star_gtf} | sort | uniq | awk '{print \$1}' | uniq -c | awk '\$1 > 1 {print \$2}' > genes_multiple_strands.txt

    grep -Fv -f transcripts_multiple_chromosomes.txt ${params.ref_star_gtf} | \\
        grep -Fv -f transcripts_multiple_strands.txt | \\
        grep -Fv -f genes_multiple_chromosomes.txt | \\
        grep -Fv -f genes_multiple_strands.txt > filtered_annotation.gtf

    # Export required vars
    export ALEVIN_FRY_HOME=.

    # simpleaf configuration
    simpleaf set-paths

    # Build reference index
    simpleaf index \\
        --fasta ${params.ref_fasta} \\
        --gtf filtered_annotation.gtf \\
        --output .
    """
}
