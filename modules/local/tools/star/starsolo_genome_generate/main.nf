process STARSOLO_INDEX {
    publishDir "${params.outdir}/genome/star_index_${ref_gtf.simpleName}_${meta.id}", mode: 'copy'
    label 'process_single_mem2'

    // Memory tracks the size of the reference being indexed, overrides process_single_mem2's flat assignments. 
    // Coefficients live in params.dynamic_memory; remove the entry to fall back to the plain label.
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.STARSOLO_INDEX, [ref_fasta, ref_gtf], task.attempt, 64) }

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/htslib_samtools_star_gawk:f196f82abbbc8871"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path ref_gtf
    path ref_fasta

    output:
    path("*"),          emit: index
    path "versions.yml", emit: versions

    script:
    // Retrieve settings from custom parameters if set, otherwise from conf/seqtech_parameters.config
    def star_genomeSAindexNbases = params.star_genomeSAindexNbases ?: params.seqtech_parameters[params.protocol].star_genomeSAindexNbases
    def star_genomeSAsparseD = params.star_genomeSAsparseD ?: params.seqtech_parameters[params.protocol].star_genomeSAsparseD

    // Derive it from task.memory * 0.85, and allow for a user override via params.star_limitGenomeGenerateRAM. 
    def genomegen_ram = params.star_limitGenomeGenerateRAM
        ?: (long) ((task.memory ? task.memory.toBytes() : 32000000000L) * 0.85)

    """
    echo "\n\n==================  GENOME INDEX STARSOLO =================="
    echo "Creating star index using GTF file: ${ref_gtf}"
    echo "--genomeSAindexNbases = ${star_genomeSAindexNbases}"
    echo "--genomeSAsparseD = ${star_genomeSAsparseD}"
    echo "--limitGenomeGenerateRAM = ${genomegen_ram} (allocation: ${task.memory})"

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat ${fastq_cDNA} | awk 'NR==2 {print length(\$0)-1; exit}' || echo "")

    echo "Generating genome index with STAR"
    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${ref_fasta} \\
        --sjdbGTFfile ${ref_gtf} \\
        --sjdbOverhang "\${sjdb_overhang}" \\
        --genomeSAsparseD ${star_genomeSAsparseD} \\
        --genomeSAindexNbases ${star_genomeSAindexNbases} \\
        --limitGenomeGenerateRAM ${genomegen_ram}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSIONS
    """
}
