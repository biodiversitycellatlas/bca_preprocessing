process MAPPING_STATS {
    publishDir "${params.outdir}/summary_results", mode: 'copy'
    label 'process_single'


    conda "${moduleDir}/environment.yml"

    input:
    val(trigger)

    output:
    path("*")

    script:
    """
    # Produce summary (.tsv) of mapping statistics
    dashboard_mappingstats.py ${params.outdir}

    # Create UMI distribution and Cell + Gene count plots
    plot_umidist_cellgenecount.R ${params.outdir}
    """
}
