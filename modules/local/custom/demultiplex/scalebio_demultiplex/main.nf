process SCALEBIO_DEMUX {
    publishDir path: { "${params.outdir}/demultiplex/scalebio/${library.id}" }, mode: 'copy',
        saveAs: { filename -> filename.endsWith('.metrics.json') || filename.contains('/unassigned.') ? filename.tokenize('/')[-1] : null }
    tag "${chunk}"

    // One lane (FASTQs) or one Ultima CRAM of one library: Python streams the reads, pigz
    // and samtools (de)compress alongside. Time follows the reads (params.dynamic_time; a
    // CRAM is counted at its FASTQ size by BcaResources), never below 1 h.
    cpus   3
    memory { 2.GB * task.attempt }
    time   { BcaResources.scaledTime(
        params.dynamic_time?.SCALEBIO_DEMUX, [fastq_cDNA, fastq_BC, fastq_I1, fastq_I2], task.attempt, 1) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pigz_samtools:016fffe7e7bd7a22' :
        'community.wave.seqera.io/library/python_pigz_samtools:17d3e4f9b8b7e970' }"

    input:
    // library: [id, ultima, samples: [[id, rt, pcr], ...]]; chunk: '<library>_<lane>', unique over
    // the run. For Ultima CRAMs the CRAM takes the cDNA slot and the other reads are empty.
    tuple val(library), val(chunk), path(fastq_cDNA, stageAs: 'in/cdna/*'), path(fastq_BC, stageAs: 'in/bc/*'), path(fastq_I1, stageAs: 'in/i1/*'), path(fastq_I2, stageAs: 'in/i2/*')
    path(barcodes)
    path(reference, stageAs: 'reference/*')     // [] unless the CRAMs are reference-compressed

    output:
    tuple val(library), val(chunk), path("out/*_cDNA.fastq.gz"), path("out/*_BC_UMI.fastq.gz"), path("out/${chunk}.metrics.json"), emit: reads
    path "versions.yml", emit: versions

    script:
    def cfg            = params.seqtech_parameters[params.protocol]
    def max_mismatches = params.scalebio_max_mismatches != null ? params.scalebio_max_mismatches : cfg.scalebio_max_mismatches
    def unassigned_arg = params.scalebio_write_unassigned ? '--write-unassigned' : ''
    def input_args
    if (library.ultima) {
        def tags = params.scalebio_ultima_tags ?: cfg.scalebio_ultima_tags
        def umi  = params.scalebio_ultima_umi ?: cfg.scalebio_ultima_umi
        input_args = "--cram ${fastq_cDNA} --ultima-tags '${tags}' --ultima-umi '${umi}' --threads ${task.cpus - 1}" +
            (reference ? " --reference ${reference}" : '')
    } else {
        input_args = "--cdna ${fastq_cDNA} --bc ${fastq_BC}" +
            (fastq_I1 ? " --i1 ${fastq_I1}" : '') + (fastq_I2 ? " --i2 ${fastq_I2}" : '')
    }
    // One line per sample, indented like the script so Nextflow's indent stripping keeps the heredoc intact
    def samples_tsv    = library.samples.collect { s -> "${s.id}\t${s.rt ?: ''}\t${s.pcr ?: ''}" }.join('\n    ')
    """
    echo -e "\\n\\n==================  ScaleBio: barcode correction and sample split  =================="
    echo "Library ${library.id}, chunk ${chunk}: samples ${library.samples.collect { s -> s.id }.join(', ')}"

    cat > samples.tsv <<'END_SAMPLES'
    ${samples_tsv}
    END_SAMPLES

    scalebio_demux.py demux \\
        --barcodes ${barcodes} \\
        --samples samples.tsv \\
        --chunk ${chunk} \\
        ${input_args} \\
        --max-mismatches ${max_mismatches} \\
        --i5-revcomp ${params.scalebio_i5_revcomp} \\
        ${unassigned_arg} \\
        --outdir out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
