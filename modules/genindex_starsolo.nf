/* GENOME INDEX STARSOLO
 * Creates the Genome Indeces required for mapping.    
*/                  

process GENINDEX_STARSOLO {   
    output:
    path("*")     
       
    script:
    def gff_arg   = task.ext.args ?: ''

    """
    echo "\n\n==================  GENOME INDEX STARSOLO =================="
    # Retrieve the first accession number
    first_fastq=\$(ls "${params.resDir}/fastq/" | head -n1)  
    echo "\${first_fastq}"

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat ${params.resDir}/fastq/\${first_fastq} 2>/dev/null | awk 'NR==2 {print length(\$0)-1; exit}' || echo "") 

    # Read configuration file
    config_file=\$(cat ${params.star_config_mkref})

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

