process STARSOLO_INDEX {
    publishDir "${params.outdir}/genome/star_index_${ref_gtf.simpleName}_${meta.id}", mode: 'copy'
    label 'process_high_memory'

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path ref_gtf
    path ref_fasta

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/htslib_samtools_star_gawk:f196f82abbbc8871"

    output:
    path("*")

    script:
    // Retrieve settings from custom parameters if set, otherwise from conf/seqtech_parameters.config
    def star_genomeSAindexNbases = params.star_genomeSAindexNbases ?: params.seqtech_parameters[params.protocol].star_genomeSAindexNbases
    def star_genomeSAsparseD = params.star_genomeSAsparseD ?: params.seqtech_parameters[params.protocol].star_genomeSAsparseD

    """
    echo "\n\n==================  GENOME INDEX STARSOLO =================="
    echo "Creating star index using GTF file: ${ref_gtf}"
    echo "--genomeSAindexNbases = ${star_genomeSAindexNbases}"
    echo "--genomeSAsparseD = ${star_genomeSAsparseD}"

    if [[ ${params.scalerna_inputFormat} == "crams" ]]; then
        # For CRAM input, use STARsolo default
        sjdb_overhang=100
    else
        # Calculate SJDB overhang using the first read from the first fastq file
        sjdb_overhang=\$(zcat ${fastq_cDNA} | awk 'NR==2 {print length(\$0)-1; exit}' || echo "")
    fi

    echo "Generating genome index with STAR"
    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${ref_fasta} \\
        --sjdbGTFfile ${ref_gtf} \\
        --sjdbOverhang "\${sjdb_overhang}" \\
        --genomeSAsparseD ${star_genomeSAsparseD} \\
        --genomeSAindexNbases ${star_genomeSAindexNbases} \\
        --limitGenomeGenerateRAM ${params.star_limitGenomeGenerateRAM}
    """
}
