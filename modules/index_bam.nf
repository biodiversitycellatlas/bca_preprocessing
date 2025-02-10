// =================  INDEX BAM FILES  ================== \\ 

process INDEX_BAM {
    publishDir "${params.resDir}/mapping_STARsolo/mapping_STARsolo_${config_name}/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    output:
    path("*.bai")

    script:
    """
    echo "\n\n==================  INDEX BAM FILES  =================="
    echo "Processing files: ${mapping_files}"

    bam_file=\$(ls ${sample_id}_${config_name}_Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"

    samtools index \${bam_file}
    """
}
