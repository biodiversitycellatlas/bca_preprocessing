process SAMTOOLS_VIEW_UNMAPPED {
    tag "${meta.id}"

    cpus   2
    memory { 1.GB * task.attempt }
    time   { 45.m * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f' :
        'community.wave.seqera.io/library/samtools:1.22.1--ab498d6f62305900' }"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.id}_unmapped.fasta"),     emit : filtered_unmapped_fasta
    path "versions.yml",                                    emit : versions

    script:
    """
    echo -e "\\n\\n==================  SAMTOOLS VIEW UNMAPPED =================="
    echo "Metadata: ${meta}"
    echo "BAM file: ${bam_file}"

    # Extract unmapped reads (include only unmapped with -f 4). samtools fasta filters
    # by itself, so the reads are not re-compressed into a BAM only to be read back
    samtools fasta -@ ${task.cpus} -f 0x4 ${bam_file} > ${meta.id}_unmapped.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
    END_VERSIONS
    """
}
