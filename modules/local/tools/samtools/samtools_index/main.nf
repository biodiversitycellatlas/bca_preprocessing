process SAMTOOLS_INDEX {
    publishDir path: { "${params.outdir}/mapping_STARsolo/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"
    label 'process_low2'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f' :
        'community.wave.seqera.io/library/samtools:1.22.1--ab498d6f62305900' }"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("*.bai"), emit: bam_index
    path "versions.yml",            emit: versions

    script:
    """
    echo -e "\\n\\n==================  SAMTOOLS INDEX  =================="
    echo "Sample ID: ${meta.id}"
    echo "Indexing BAM file: ${bam_file}"

    samtools index -@ ${task.cpus} ${bam_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
    END_VERSIONS
    """
}
