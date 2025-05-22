/* GENOME INDEX STARSOLO
 * Creates the Genome Indeces required for mapping.    
*/                  

process GENINDEX_STARSOLO {
    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)
    path ref_gtf

    output:
    path("*")  
       
    script:
    def gff_arg   = task.ext.args ?: ''

    // Retrieve star_index settings from conf/seqtech_parameters.config
    def star_index_settings = params.seqtech_parameters[params.protocol].star_index
    def star_index_args = star_index_settings instanceof List ? star_index_settings.join(' ') : star_index_settings

    """
    echo "\n\n==================  GENOME INDEX STARSOLO =================="
    echo "Creating star index using GTF file: ${ref_gtf}" 

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat ${fastq_cDNA} | awk 'NR==2 {print length(\$0)-1; exit}' || echo "") 

    echo "Generating genome index with STAR"
    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${params.ref_fasta} \\
        --sjdbGTFfile ${ref_gtf} \\
        --sjdbOverhang "\${sjdb_overhang}" \\
        ${gff_arg} \\
        ${star_index_args}
    """
}

