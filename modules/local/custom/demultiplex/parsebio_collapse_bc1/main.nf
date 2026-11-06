process PARSEBIO_COLLAPSE_BC1 {
    publishDir path: { "${params.outdir}/demultiplex/bc1_collapse/${meta.id}" }, mode: 'copy', pattern: '*_bc1_collapse.log'
    tag "${meta.id}"

    // One pass over the barcode read; time follows its size.
    cpus   1
    memory { 2.GB * task.attempt }
    time   { BcaResources.scaledTime(
        params.dynamic_time?.PARSEBIO_COLLAPSE_BC1, [fastq_BC_UMI], task.attempt, 4) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930' :
        'community.wave.seqera.io/library/python:3.14.2--0562a0df3245213a' }"

    input:
    tuple val(meta), path(fastq_BC_UMI)
    path(whitelist, stageAs: 'splitwells_whitelist/*')

    output:
    tuple val(meta), path("*_R2.bc1collapsed.fastq.gz"), emit: reads
    path "parse_bc1_T_whitelist.txt",                    emit: bc1_whitelist
    path "*_bc1_collapse.log",                           emit: log
    path "versions.yml",                                 emit: versions

    script:
    """
    echo -e "\\n\\n==================  Parse Biosciences: merge round-1 barcodes per well  =================="
    echo "Processing sample: ${meta.id}"

    # Rewrite bc1 (R2 positions 51-58) to the well's oligo-dT barcode
    parsebio_collapse_bc1.py \\
        --csv ${whitelist} \\
        --fq2 ${fastq_BC_UMI} \\
        --sample_id ${meta.id} \\
        --output . \\
        --barcode_start 50 \\
        --barcode_end 58 \\
        --max_edit_dist 2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
