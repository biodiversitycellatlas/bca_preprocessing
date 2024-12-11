// =============  GENOME INDEX STARSOLO  ============= \\ 
// Creates the Genome Indeces required for mapping.    \\
// Two folders will be created within /genome/,        \\
// genome_index: contains genome index for STARsolo    \\
// parse_refgenome: see next process                   \\

process GENINDEX_STARSOLO {   
    publishDir "${params.resDir}/genome/genome_index_${config_name}", mode: 'symlink'
    debug true
    label 'big_mem'

    input:
    path(ref_star_gff)
    path(config_mkref)
    val config_name

    output:
    path("*")     
       
    script:
    """
    echo "\n\n==================  GENOME INDEX STARSOLO ${config_name} =================="
    # Retrieve the first accession number
    first_accs=\$(head -1 ${params.accessions})
    first_fastq="${params.resDir}/fastq/\${first_accs}*1*.fastq.gz"  

    echo "\${first_accs}"
    echo "\${first_fastq}"

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat \${first_fastq} 2>/dev/null | awk 'NR==2 {print length(\$0)-1; exit}' || echo "") 

    # Read configuration file
    config_file=\$(cat ${config_mkref})

    echo "Generating genome index with STAR"
    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${params.ref_fasta} \\
        --sjdbGTFfile ${ref_star_gff} \\
        --sjdbOverhang "\${sjdb_overhang}" \\
        --genomeSAindexNbases 12 \\
        \${config_file}

    # Only add this parameter for GFF files
    # --sjdbGTFtagExonParentTranscript transcript_id \\
    """
}

