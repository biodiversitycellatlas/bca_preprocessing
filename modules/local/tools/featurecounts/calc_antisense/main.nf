process CALC_ANTISENSE {
    publishDir "${params.outdir}/rRNA_mtDNA", mode: 'copy'
    tag "${meta.id}"
    label 'process_single_long2'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/samtools_subread:f5fd17c543add0fd"

    input:
    tuple val(meta), path(bam_file)
    file(ref_gtf)

    output:
    path("${meta.id}_antisense_metrics.txt"), emit: antisense_metrics
    path "versions.yml",                      emit: versions

    script:
    """
    echo "\n\n==================  ANTISENSE READS =================="
    echo "Sample ID: ${meta.id}"
    echo "BAM file: ${bam_file}"
    echo "GTF: ${ref_gtf}"
    echo "STARsolo strand: ${params.star_soloStrand}"

    calculate_antisense.sh \\
        ${bam_file} \\
        ${ref_gtf} \\
        ${params.star_soloStrand} \\
        ${meta.id}_antisense_metrics.txt \\
        ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | sed 's/featureCounts v//')
    END_VERSIONS
    """
}
