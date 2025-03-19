process INDEX_BAM {
    publishDir "${params.resDir}/mapping_STARsolo/mapping_STARsolo/${sample_id}", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(mapping_files)

    output:
    path("${sample_id}_Aligned.sortedByCoord.out.bam.bai")

    script:
    """
    echo "\n\n==================  INDEX BAM FILES  =================="
    echo "Sample ID: ${sample_id}"
    echo "Processing files: ${mapping_files}"

    bam_file=\$(ls ${sample_id}_Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"

    samtools index \${bam_file}
    """
}
