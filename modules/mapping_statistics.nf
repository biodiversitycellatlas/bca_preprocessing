process MAPPING_STATS {
    publishDir "${params.output_dir}/summary_results", mode: 'copy', overwrite: true
    tag "${sample_ids}"
    debug true

    input:
    val(trigger)
    val(sample_ids)

    output:
    path("*")

    script:
    """
    echo "Results directory: ${params.output_dir}"
    echo "Sample IDs: ${sample_ids}"

    # Produce summary (.tsv) of mapping statistics
    sbatch ${params.code_dir}/bin/mapping_statistics.sh ${params.output_dir} ${params.ref_gtf}

    # Create UMI distribution and Cell + Gene count plots
    Rscript ${params.code_dir}/bin/plot_umidist_cellgenecount.R ${params.output_dir} ${sample_ids}
    """
}
