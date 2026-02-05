process KRONA {
    publishDir "${params.outdir}/summary_results", mode: 'copy'
    label 'process_single'


    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/krona:2.8--8fc5aef4acd456a6"

    input:
    val(trigger)

    output:
    path("*")

    script:
    """
    # Installing/updating krona database
    ktUpdateTaxonomy.sh

    # Running Krona on Kraken reports
    # -t is set to 7 as kraken reports are created with --report-minimizer-data
    ktImportTaxonomy -t 7 -m 3 -o multi-krona.html ${params.outdir}/*/*_taxonomy.txt
    """
}
