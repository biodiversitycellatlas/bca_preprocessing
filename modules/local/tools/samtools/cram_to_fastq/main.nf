process CRAM_TO_FASTQ {
    tag "${meta.id}"

    // samtools decodes and compresses with its threads. Time follows the CRAM
    // (params.dynamic_time), which BcaResources counts at its FASTQ size, never below 1 h.
    cpus   4
    memory { 2.GB * task.attempt }
    time   { BcaResources.scaledTime(
        params.dynamic_time?.CRAM_TO_FASTQ, [cram], task.attempt, 1) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pigz_samtools:016fffe7e7bd7a22' :
        'community.wave.seqera.io/library/python_pigz_samtools:17d3e4f9b8b7e970' }"

    input:
    // One unaligned CRAM/BAM; meta: [id: its name, path: its full path, to give the FASTQs
    // back to every sample that reads it]
    tuple val(meta), path(cram)
    path(reference, stageAs: 'reference/*')     // [] unless the CRAMs are reference-compressed

    output:
    tuple val(meta), path("*_R1_001.fastq.gz"), path("*_R2_001.fastq.gz"), path("*_I{1,2}_001.fastq.gz", arity: '0..2'), emit: fastqs
    path "versions.yml", emit: versions

    script:
    // Named like bcl-convert output, so the lanes and index reads are recognised downstream
    def stem = cram.name.replaceFirst(/\.(cram|bam|sam)$/, '')
    def ref_arg = reference ? "--reference ${reference}" : ''
    def index_format = params.cram_index_format ?: 'auto'
    """
    echo -e "\\n\\n==================  CRAM to FASTQ: ${cram.name}  =================="

    # The index reads are in the BC (sequence) and QT (quality) tags, the two indices of a
    # dual-indexed read separated by a non-base character (e.g. 'ACGT-TTGG'). 'auto' counts
    # them in the first 10,000 reads; tags without a separator need 'cram_index_format'.
    index_format="${index_format}"
    if [ "\$index_format" = "auto" ]; then
        n_index=\$(python3 - ${cram} ${reference ?: ''} <<'PY'
    import re, subprocess, sys
    cmd = ["samtools", "view"] + (["-T", sys.argv[2]] if len(sys.argv) > 2 else []) + [sys.argv[1]]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    n = 0
    for i, line in enumerate(p.stdout):
        for t in line.rstrip("\\n").split("\\t")[11:]:
            if t.startswith("BC:Z:"):
                n = max(n, len([x for x in re.split(r"[^ACGTNacgtn]+", t[5:]) if x]))
        if i >= 9999:
            break
    p.kill()
    print(n)
    PY
        )
        case "\$n_index" in
            0) index_format="" ;;
            1) index_format="i*" ;;
            *) index_format="i*i*" ;;
        esac
    fi
    echo "Index reads: \${index_format:-none}"

    index_args=()
    if [ -n "\$index_format" ]; then
        index_args=(--barcode-tag BC --quality-tag QT --index-format "\$index_format" --i1 ${stem}_I1_001.fastq.gz)
        only_i="\${index_format//[^i]/}"
        if [ "\${#only_i}" -ge 2 ]; then
            index_args+=(--i2 ${stem}_I2_001.fastq.gz)
        fi
    fi

    # Reads that are not one of a pair (single-end, e.g. an Ultima CRAM) would be lost
    samtools fastq -@ ${task.cpus} ${ref_arg} -n "\${index_args[@]}" \\
        -1 ${stem}_R1_001.fastq.gz -2 ${stem}_R2_001.fastq.gz \\
        -0 unpaired_other.fastq -s unpaired_singletons.fastq \\
        ${cram}
    if [ -s unpaired_other.fastq ] || [ -s unpaired_singletons.fastq ]; then
        echo "ERROR: ${cram.name} holds reads that are not one of a read pair (single-end data?)." >&2
        echo "Illumina-style CRAMs need both reads of every pair; Ultima CRAMs take 'cram_type = \\"ultima\\"' (ScaleBio)." >&2
        exit 1
    fi
    rm -f unpaired_other.fastq unpaired_singletons.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
