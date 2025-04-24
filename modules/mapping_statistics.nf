process MAPPING_STATS {
    publishDir "${params.output_dir}/summary_results", mode: 'copy', overwrite: true
    tag "${meta.id}"
    debug true

    input:
    val(trigger)
    val(meta)

    output:
    path("*")

    script:
    """
    echo "Results directory: ${params.output_dir}"
    echo "Sample IDs: ${meta}"

    # Produce summary (.tsv) of mapping statistics
    sbatch ${launchDir}/bin/mapping_statistics.sh ${params.output_dir} ${params.ref_gtf}

    # Create UMI distribution and Cell + Gene count plots
    Rscript ${launchDir}/bin/plot_umidist_cellgenecount.R ${params.output_dir} ${meta.id}
    """
}
