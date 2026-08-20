process SAMTOOLS_MERGE {
    tag "merging_bams"
    label 'process_single_mem2'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"

    input:
    path(bams)

    output:
    path "merged_genome.bam" ,      emit: merged_bam
    path "merged_genome.bam.bai" ,  emit: merged_bai
    path "versions.yml",            emit: versions

    script:
    """
    input_count=\$(echo ${bams} | wc -w)

    if [ "\$input_count" -eq 1 ]; then
        echo "Only one BAM found. Skipping merge and symlinking..."
        ln -s ${bams} merged_genome.bam
    else
        echo "Merging \$input_count BAM files..."
        samtools merge -@ ${task.cpus} merged_genome.bam ${bams}
    fi

    samtools index -@ ${task.cpus} merged_genome.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
    END_VERSIONS
    """
}
