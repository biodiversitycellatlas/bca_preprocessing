process CALC_MT_RRNA {
    publishDir "${params.outdir}/rRNA_mtDNA", mode: 'copy'
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/samtools_subread:f5fd17c543add0fd"

    input:
    tuple val(meta), path(bam_file)
    file(bam_index)
    file(ref_gtf)
    file(rrna_gtf)

    output:
    path("${meta.id}_mt_rrna_metrics.txt"), emit: mt_rrna_metrics
    path "versions.yml",                    emit: versions

    script:
    def rrna_arg = rrna_gtf ? "--rrna-gtf ${rrna_gtf} " : ""
    """
    echo "\n\n==================  CALCULATION rRNA & mtDNA =================="
    echo "Sample ID: ${meta.id}"
    echo "BAM file: ${bam_file}"
    echo "GTF: ${ref_gtf}"
    echo "Added rRNA reference: ${rrna_gtf ?: 'none'}"

    calculate_rrna_mtdna.sh \\
        ${rrna_arg}\\
        ${bam_file} \\
        ${meta.id}_mt_rrna_metrics.txt \\
        ${ref_gtf} \\
        "${params.mt_contig}"

    echo "Created file: ${meta.id}_mt_rrna_metrics.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
        subread: \$(featureCounts -v 2>&1 | sed 's/featureCounts v//')
    END_VERSIONS
    """
}
