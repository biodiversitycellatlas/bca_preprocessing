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

    # If mapping_software is STARsolo or both, create UMI distribution and Cell + Gene count plots
    if [ "${params.mapping_software}" == "starsolo" ] || [ "${params.mapping_software}" == "both" ]; then
        plot_umidist_cellgenecount.R ${params.outdir}
    fi
    """
}
