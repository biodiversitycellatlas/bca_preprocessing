process INDEX_BAM {
    publishDir "${params.output_dir}/mapping_STARsolo/${meta.id}", mode: 'copy'
    tag "${meta.id}"

    input:
    tuple val(meta), path(mapping_files)

    output:
    path("${meta.id}_Aligned.sortedByCoord.out.bam.bai")

    script:
    """
    echo "\n\n==================  INDEX BAM FILES  =================="
    echo "Sample ID: ${meta}"
    echo "Processing files: ${mapping_files}"

    bam_file=\$(ls ${meta.id}_Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"

    samtools index \${bam_file}
    """
}
