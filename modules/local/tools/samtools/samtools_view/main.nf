process SAMTOOLS_VIEW {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"

    input:
    tuple val(meta), path(mapping_files)

    output:
    tuple val(meta), path("${meta.id}_Aligned.filtered.sorted.bam"),     emit : filtered_bam
    tuple val(meta), path("${meta.id}_Aligned.filtered.sorted.bam.bai"), emit : filtered_bam_index
    tuple val(meta), path("mapreads.txt"),                               emit : mapreads

    script:
    """
    echo "\n\n==================  SAMTOOLS VIEW  =================="
    echo "Sample ID: ${meta}"
    echo "Processing files: ${mapping_files}"

    # Remove unmapped reads from the BAM file
    samtools view -b -F 4 ${meta.id}_Aligned.sortedByCoord.out.bam | samtools sort -o ${meta.id}_Aligned.filtered.sorted.bam
    samtools index ${meta.id}_Aligned.filtered.sorted.bam

    # Calculate mapped reads
    samtools view -F 260 ${meta.id}_Aligned.filtered.sorted.bam | wc -l > mapreads.txt
    """
}
