process SAMTOOLS_INDEX {
    publishDir "${params.outdir}/mapping_STARsolo/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("*.bai"), emit: bam_index

    script:
    """
    echo "\n\n==================  SAMTOOLS INDEX  =================="
    echo "Sample ID: ${meta.id}"
    echo "Indexing BAM file: ${bam_file}"

    samtools index -@ ${task.cpus} ${bam_file}
    """
}
