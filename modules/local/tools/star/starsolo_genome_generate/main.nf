process STARSOLO_INDEX {
    publishDir path: { "${params.outdir}/genome/star_index_${ref_gtf.simpleName}${index_suffix}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_single_mem2'

    // Memory tracks the size of the reference being indexed, overrides process_single_mem2's flat assignments. 
    // Coefficients live in params.dynamic_memory; remove the entry to fall back to the plain label.
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.STARSOLO_INDEX, [ref_fasta, ref_gtf], task.attempt, 64) }

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/htslib_samtools_star_gawk:f196f82abbbc8871"

    // One index per run: fastq_cDNA holds the cDNA reads of every sample, only used for
    // the read length. Staged into numbered dirs so equal filenames across samples don't clash.
    input:
    path(fastq_cDNA, stageAs: 'cdna_?/*')
    path ref_gtf
    path ref_fasta
    val index_suffix

    output:
    path("GenomeDir"),  emit: index
    path "versions.yml", emit: versions

    script:
    // Retrieve settings from custom parameters if set, otherwise from conf/seqtech_parameters.config
    def star_genomeSAindexNbases = params.star_genomeSAindexNbases ?: params.seqtech_parameters[params.protocol].star_genomeSAindexNbases
    def star_genomeSAsparseD = params.star_genomeSAsparseD ?: params.seqtech_parameters[params.protocol].star_genomeSAsparseD

    // Derive it from task.memory * 0.85, and allow for a user override via params.star_limitGenomeGenerateRAM. 
    def genomegen_ram = params.star_limitGenomeGenerateRAM
        ?: ((task.memory ? task.memory.toBytes() : 32000000000L) * 0.85) as long

    """
    echo -e "\\n\\n==================  GENOME INDEX STARSOLO =================="
    echo "Creating star index using GTF file: ${ref_gtf}"
    echo "--genomeSAindexNbases = ${star_genomeSAindexNbases}"
    echo "--genomeSAsparseD = ${star_genomeSAsparseD}"
    echo "--limitGenomeGenerateRAM = ${genomegen_ram} (allocation: ${task.memory})"

    # SJDB overhang = max read length - 1, taken from the first read of every sample's cDNA fastq
    sjdb_overhang=""
    for fq in cdna_*/*; do
        len=\$(zcat "\$fq" | awk 'NR==2 {print length(\$0)-1; exit}' || echo "")
        if [ -n "\$len" ] && { [ -z "\$sjdb_overhang" ] || [ "\$len" -gt "\$sjdb_overhang" ]; }; then
            sjdb_overhang=\$len
        fi
    done
    echo "--sjdbOverhang = \${sjdb_overhang}"

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
