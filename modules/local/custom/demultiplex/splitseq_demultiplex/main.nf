process SPLITSEQ_DEMUX {
    publishDir path: { "${params.outdir}/demultiplex/demux_splitseq/${meta.id}" }, mode: 'copy', pattern: '*_demux.log'
    tag "${meta.id}"

    // One pass over both reads; time follows their size (params.dynamic_time)
    cpus   1
    memory { 2.GB * task.attempt }
    time   { BcaResources.scaledTime(
        params.dynamic_time?.SPLITSEQ_DEMUX, [fastq_cDNA, fastq_BC_UMI], task.attempt, 4) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930' :
        'community.wave.seqera.io/library/python:3.14.2--0562a0df3245213a' }"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path(wells, stageAs: 'splitseq_wells/*')
    val(bc1_start)                              // round-1 barcode in read 2, 0-based start
    val(bc1_end)                                // and exclusive end

    output:
    tuple val(meta), path("*group*_R1*"), path("*group*_R2*"), path(fastq_indices), path(input_file), emit: splitted_files
    path "*_demux.log",  emit: log
    path "versions.yml", emit: versions

    script:
    """
    echo -e "\\n\\n==================  SPLiT-seq: split the reads on the RT wells  =================="
    echo "Sample: ${meta.id}, wells (rt): ${meta.rt}"

    # Keeps the pairs whose round-1 barcode (smallest Levenshtein distance, at most 2)
    # lies in one of the sample's wells, oligo-dT and random-hexamer barcodes alike
    parsebio_custom_demux.py \\
        --sample_id ${meta.id} \\
        --fq1 ${fastq_cDNA} \\
        --fq2 ${fastq_BC_UMI} \\
        --whitelist ${wells} \\
        --group ${meta.id} '${meta.rt}' \\
        --output . \\
        --barcode_start ${bc1_start} \\
        --barcode_end ${bc1_end} \\
        --max_edit_dist 2 \\
        > ${meta.id}_demux.log
    cat ${meta.id}_demux.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
