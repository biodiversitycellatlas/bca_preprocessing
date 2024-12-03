// =============  GENOME INDEX STARSOLO  ============= \\ 
// Creates the Genome Indeces required for mapping.    \\
// Two folders will be created within /genome/,        \\
// genome_index: contains genome index for STARsolo    \\
// parse_refgenome: see next process                   \\

process GENINDEX_STARSOLO {   
    publishDir "${params.resDir}/genome/genome_index", mode: 'symlink'
    debug true

    input:
    file(ref_star_gtf)

    output:
    path("*")     
       
    script:
    """
    echo "\n\n==================  GENOME INDEX STARSOLO  =================="
    # Retrieve the first accession number
    first_accs=\$(head -1 ${params.dataDir}/accession_lists/${params.species}_${params.seqTech}_accessions.txt)
    first_fastq="${params.resDir}/fastq/\${first_accs}*1*.fastq.gz"  

    echo "\${first_accs}"
    echo "\${first_fastq}"

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat \${first_fastq} 2>/dev/null | awk 'NR==2 {print length(\$0)-1; exit}' || echo "") 

    echo "\${sjdb_overhang}"

    echo "Generating genome index with STAR"
    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${params.ref_fasta} \\
        --sjdbGTFfile ${ref_star_gtf} \\
        --sjdbOverhang "\${sjdb_overhang}" \\
        --genomeSAindexNbases 12
    """
}

