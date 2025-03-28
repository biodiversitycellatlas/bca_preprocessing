process MAPPING_STATS {
    publishDir "${params.resDir}/summary_results", mode: 'copy', overwrite: true
    tag "${sample_ids}"
    debug true

    input:
    val(trigger)
    val(sample_ids)

    output:
    path("*")

    script:
    """
    echo "Results directory: ${params.resDir}"
    echo "Sample IDs: ${sample_ids}"

    # Produce summary (.tsv) of mapping statistics
    sbatch ${params.baseDir}/bin/mapping_statistics.sh ${params.resDir}

    # Create UMI distribution and Cell + Gene count plots
    Rscript ${params.baseDir}/bin/plot_umidist_cellgenecount.R ${params.resDir} ${sample_ids}
    """
}
