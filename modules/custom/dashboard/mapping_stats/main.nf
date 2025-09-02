process MAPPING_STATS {
    publishDir "${params.outdir}/summary_results", mode: 'copy'
    label 'process_single'
    debug true

    conda "${moduleDir}/environment.yml"

    input:
    val(trigger)

    output:
    path("*")

    script:
    """
    # Produce summary (.tsv) of mapping statistics
    python3 ${moduleDir}/bin/dashboard_mappingstats.py ${params.outdir}

    # Create UMI distribution and Cell + Gene count plots
    Rscript ${moduleDir}/bin/plot_umidist_cellgenecount.R ${params.outdir}
    """
}
