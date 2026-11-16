process STARSOLO_INDEX {
    publishDir path: { "${params.outdir}/genome/star_index_${ref_gtf.simpleName}${index_suffix}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }

    // Memory and time track the size of the reference being indexed. CPUs step up from 2
    // above a reference size in params.dynamic_cpus: a small genome indexes on about one core.
    cpus   { BcaResources.tieredCpus(
        params.dynamic_cpus?.STARSOLO_INDEX, [ref_fasta, ref_gtf], task.attempt, 2) }
    time   { BcaResources.scaledTime(
        params.dynamic_time?.STARSOLO_INDEX, [ref_fasta, ref_gtf], task.attempt, 2) }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.STARSOLO_INDEX, [ref_fasta, ref_gtf], task.attempt, 10) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/htslib_samtools_star_gawk:f196f82abbbc8871' :
        'community.wave.seqera.io/library/htslib_samtools_star_gawk:782fa39fdeb1558c' }"

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

    // A user override of --genomeChrBinNbits; empty = derive it from the reference below
    def chr_bin_override = params.star_genomeChrBinNbits ?: ''

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

    # --genomeChrBinNbits = min(18, floor(log2(max(genome length / sequences, read length)))),
    # the STAR manual's scaling: every sequence takes a whole number of bins, so a genome with
    # many short contigs wastes memory on the default of 18
    if [ -n "${chr_bin_override}" ]; then
        chr_bin_nbits=${chr_bin_override}
    else
        read genome_len n_refs < <(awk '/^>/ {n++; next} {sub(/\\r\$/, ""); len += length(\$0)} END {print len+0, n+0}' ${ref_fasta})
        read_len=\$(( \${sjdb_overhang:-0} + 1 ))
        chr_bin_nbits=\$(awk -v g="\$genome_len" -v n="\$n_refs" -v r="\$read_len" 'BEGIN {
            if (n < 1) { print 18; exit }
            x = g / n; if (r > x) x = r
            b = int(log(x) / log(2)); if (b > 18) b = 18; if (b < 1) b = 1
            print b }')
        echo "Reference: \${genome_len} bp in \${n_refs} sequences, read length \${read_len}"
    fi
    echo "--genomeChrBinNbits = \${chr_bin_nbits}"

    echo "Generating genome index with STAR"
    STAR --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeFastaFiles ${ref_fasta} \\
        --sjdbGTFfile ${ref_gtf} \\
        --sjdbOverhang "\${sjdb_overhang}" \\
        --genomeSAsparseD ${star_genomeSAsparseD} \\
        --genomeSAindexNbases ${star_genomeSAindexNbases} \\
        --genomeChrBinNbits "\${chr_bin_nbits}" \\
        --limitGenomeGenerateRAM ${genomegen_ram}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSIONS
    """
}
