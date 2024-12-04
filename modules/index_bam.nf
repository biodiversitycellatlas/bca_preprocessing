// =================  INDEX BAM FILES  ================== \\ 

process INDEX_BAM {
    publishDir "${params.resDir}/mapping_parsepipe/${sample_id}", mode: 'symlink'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    output:
    tuple val(sample_id), val(config_name), path("*.bam.bai")

    script:
    """
    echo "\n\n==================  INDEX BAM FILES  =================="
    echo "Processing files: ${mapping_files}"

    bam_file=\$(ls *Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"

    samtools index \${bam_file}
    """
}
