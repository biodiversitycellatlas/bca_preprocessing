process SPLITSEQ_COLLAPSE_BC1 {
    publishDir path: { "${params.outdir}/demultiplex/bc1_collapse/${meta.id}" }, mode: 'copy', pattern: '*_bc1_collapse.log'
    tag "${meta.id}"

    // One pass over the barcode read; time follows its size
    cpus   1
    memory { 1.GB * task.attempt }
    time   { BcaResources.scaledTime(
        params.dynamic_time?.SPLITSEQ_COLLAPSE_BC1, [fastq_BC_UMI], task.attempt, 4) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930' :
        'community.wave.seqera.io/library/python:3.14.2--0562a0df3245213a' }"

    input:
    tuple val(meta), path(fastq_BC_UMI)
    path(wells, stageAs: 'splitseq_wells/*')
    val(bc1_start)                              // round-1 barcode in read 2, 0-based start
    val(bc1_end)                                // and exclusive end

    output:
    tuple val(meta), path("*_R2.bc1collapsed.fastq.gz"), emit: reads
    path "splitseq_bc1_T_whitelist.txt",                 emit: bc1_whitelist
    path "*_bc1_collapse.log",                           emit: log
    path "versions.yml",                                 emit: versions

    script:
    """
    echo -e "\\n\\n==================  SPLiT-seq: merge the round-1 barcodes per RT well  =================="
    echo "Processing sample: ${meta.id}"

    # Rewrite bc1 to the well's oligo-dT barcode, so a cell's oligo-dT and random-hexamer
    # reads share one cell barcode. Stops when most reads are too short to hold bc1.
    parsebio_collapse_bc1.py \\
        --csv ${wells} \\
        --fq2 ${fastq_BC_UMI} \\
        --sample_id ${meta.id} \\
        --output . \\
        --barcode_start ${bc1_start} \\
        --barcode_end ${bc1_end} \\
        --max_edit_dist 2 \\
        --whitelist_name splitseq_bc1_T_whitelist.txt \\
        --max_too_short 0.5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
