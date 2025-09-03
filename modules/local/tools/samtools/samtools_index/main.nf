process SAMTOOLS_INDEX {
    publishDir "${params.outdir}/mapping_STARsolo/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(mapping_files)

    output:
    path("${meta.id}_Aligned.sortedByCoord.out.bam.bai")

    script:
    """
    echo "\n\n==================  SAMTOOLS INDEX  =================="
    echo "Sample ID: ${meta}"
    echo "Processing files: ${mapping_files}"

    bam_file=\$(ls ${meta.id}_Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"

    samtools index \${bam_file}
    """
}
