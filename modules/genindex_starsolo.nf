// =============  GENOME INDEX STARSOLO  ============= \\ 
// Creates the Genome Indeces required for mapping.    \\
// Two folders will be created within /genome/,        \\
// genome_index: contains genome index for STARsolo    \\
// parse_refgenome: see next process                   \\

process GENINDEX_STARSOLO {   
    input:
    path(ref_star_gtf)
    path(config_mkref)
    val config_name

    output:
    path("*")     
       
    script:
    def gff_arg   = task.ext.args ?: ''          // If ext.args is defined assign it to gff_arg

    """
    echo "\n\n==================  GENOME INDEX STARSOLO ${config_name} =================="
    # Retrieve the first accession number
    first_accs=\$(head -1 ${params.accessions})
    first_fastq="${params.resDir}/fastq/\${first_accs}*R1*.fastq.gz"  

    echo "\${first_accs}"
    echo "\${first_fastq}"

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat \${first_fastq} 2>/dev/null | awk 'NR==2 {print length(\$0)-1; exit}' || echo "") 

    # Read configuration file
    config_file=\$(cat ${config_mkref})

    echo "Generating genome index with STAR"
    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${params.ref_fasta} \\
        --sjdbGTFfile ${params.ref_star_gtf} \\
        ${gff_arg} \\
        --sjdbOverhang "\${sjdb_overhang}" \\
        --genomeSAindexNbases 12 \\
        \${config_file}
    """
}

