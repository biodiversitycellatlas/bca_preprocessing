process STARSOLO_INDEX {
    publishDir "${params.outdir}/genome/star_index_${meta.id}", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path ref_gtf

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data' :
        'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4' }"

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

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat ${fastq_cDNA} | awk 'NR==2 {print length(\$0)-1; exit}' || echo "")

    echo "Generating genome index with STAR"
    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${params.ref_fasta} \\
        --sjdbGTFfile ${ref_gtf} \\
        --sjdbOverhang "\${sjdb_overhang}" \\
        --genomeSAsparseD ${star_genomeSAsparseD} \\
        --genomeSAindexNbases ${star_genomeSAindexNbases}
    """
}
