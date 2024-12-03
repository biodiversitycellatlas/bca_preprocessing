// =================  INDEX BAM FILES  ================== \\ 

process INDEX_BAM {
    publishDir "${params.resDir}/mapping_parsepipe/${sample_id}", mode: 'symlink'
    debug true

    input:
    path(mapping_files)

    output:
    file("Aligned.sortedByCoord.out.bam.bai")

    script:
    """
    echo "\n\n==================  INDEX BAM FILES  =================="
    echo "Processing files: ${mapping_files}"

    samtools index ${mapping_files}/Aligned.sortedByCoord.out.bam
    """
}
