process SAMTOOLS_VIEW_MAPPED {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.id}_Aligned.filtered.sorted.bam"),     emit : filtered_mapped_bam
    tuple val(meta), path("${meta.id}_Aligned.filtered.sorted.bam.bai"), emit : filtered_mapped_bai
    tuple val(meta), path("mapreads.txt"),                               emit : mapreads
    path "versions.yml",                                                 emit : versions

    script:
    """
    echo "\n\n==================  SAMTOOLS VIEW MAPPED =================="
    echo "Sample ID: ${meta}"
    echo "Processing files: ${bam_file}"

    # Extract mapped reads (exclude unmapped with -F 4), sort, and index
    samtools view --threads ${task.cpus} -b -F 4 ${bam_file} | samtools sort -o ${meta.id}_Aligned.filtered.sorted.bam
    samtools index ${meta.id}_Aligned.filtered.sorted.bam

    # Calculate mapped reads
    samtools view --threads ${task.cpus} -F 260 ${meta.id}_Aligned.filtered.sorted.bam | wc -l > mapreads.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
    END_VERSIONS
    """
}
