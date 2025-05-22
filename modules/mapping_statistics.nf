process MAPPING_STATS {
    publishDir "${params.output_dir}/summary_results", mode: 'copy'
    debug true

    input:
    val(trigger)

    output:
    path("*")

    script:
    """
    # MultiQC report
    multiqc ${params.output_dir} --config ${launchDir}/bin/multiqc_config.yml

    # Produce summary (.tsv) of mapping statistics
    sbatch ${launchDir}/bin/mapping_statistics.sh ${params.output_dir} ${params.ref_gtf}

    # Create UMI distribution and Cell + Gene count plots
    # Rscript ${launchDir}/bin/plot_umidist_cellgenecount.R ${params.output_dir}
    """
}
