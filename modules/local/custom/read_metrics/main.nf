process CALC_READ_METRICS {
    publishDir path: { "${params.outdir}/rRNA_mtDNA" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"
    label 'process_single_long2'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/pysam_samtools_bc_python_pruned:82a1e27e868113f0"

    input:
    tuple val(meta), path(bam_file), path(filtered_matrix_dir), path(cellreads_stats)
    path(ref_gtf)
    path(rrna_gtf)

    output:
    path("${meta.id}_mt_rrna_metrics.txt"),                   emit: mt_rrna_metrics
    path("${meta.id}_antisense_metrics.txt"),                 emit: antisense_metrics, optional: true
    tuple val(meta), path("${meta.id}_barcode_reads.tsv.gz"), emit: barcode_reads
    path "versions.yml",                                      emit: versions

    script:
    // The added features are the rRNA reference: every read on their contigs counts as rRNA
    def rrna_arg      = rrna_gtf ? "--rrna-gtf ${rrna_gtf}" : ""

    // The filtered matrix's barcodes add the called-cell rows; runs without one
    // (bam_only, or no cell calling) report those rows as N/A
    def cells_arg     = filtered_matrix_dir ? "--cell-barcodes ${filtered_matrix_dir}/barcodes.tsv" : ""
    
    // Antisense comes from STARsolo's own read flags; without CellReads.stats it is skipped
    def cellreads_arg = cellreads_stats ? "--cellreads ${cellreads_stats}" : ""
    """
    echo -e "\\n\\n==================  READ METRICS: rRNA, mtDNA & ANTISENSE =================="
    echo "Sample ID: ${meta.id}"
    echo "BAM file: ${bam_file}"
    echo "GTF: ${ref_gtf}"
    echo "Added rRNA reference: ${rrna_gtf ?: 'none'}"
    echo "Called-cell barcodes: ${filtered_matrix_dir ? "${filtered_matrix_dir}/barcodes.tsv" : 'none'}"
    echo "CellReads.stats: ${cellreads_stats ?: 'none'}"
    echo "STARsolo strand: ${params.star_soloStrand}"

    calculate_read_metrics.py \\
        --bam ${bam_file} \\
        --gtf ${ref_gtf} \\
        --mt-contig "${params.mt_contig}" \\
        --strand ${params.star_soloStrand} \\
        ${rrna_arg} \\
        ${cells_arg} \\
        ${cellreads_arg} \\
        --out-mt-rrna ${meta.id}_mt_rrna_metrics.txt \\
        --out-antisense ${meta.id}_antisense_metrics.txt \\
        --out-barcode-reads ${meta.id}_barcode_reads.tsv.gz \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """
}
