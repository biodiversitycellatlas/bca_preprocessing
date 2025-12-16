process CALC_MT_RRNA {
    publishDir "${params.outdir}/rRNA_mtDNA", mode: 'copy'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/samtools_subread:f5fd17c543add0fd"

    input:
    tuple val(meta), path(mapping_files)
    file(bam_index)

    output:
    file("${meta.id}_mt_rrna_metrics.txt")

    script:
    """
    echo "\n\n==================  CALCULATION rRNA & mtDNA =================="
    echo "Sample ID: ${meta}"

    bam_file=\$(ls ${meta.id}_Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"

    calculate_rrna_mtdna.sh \\
        \${bam_file} \\
        ${meta.id}_mt_rrna_metrics.txt \\
        ${params.ref_gtf} \\
        ${params.grep_rrna} \\
        ${params.mt_contig}

    echo "Created file: ${meta.id}_mt_rrna_metrics.txt"
    """
}
