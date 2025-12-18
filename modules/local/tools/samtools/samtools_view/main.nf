process SAMTOOLS_VIEW {
    publishDir "${params.outdir}/mapping_STARsolo/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'quay.io/biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(mapping_files)

    output:
    path("${meta.id}_Aligned.filtered.sorted.bam"),     emit : bam_file
    path("${meta.id}_Aligned.filtered.sorted.bam.bai"), emit : bam_index
    path("mapreads.txt"),                               emit : mapreads

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
