process KRONA {
    publishDir "${params.outdir}/kraken", mode: 'copy'
    label 'process_single'
    debug true

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/krona:2.8--8fc5aef4acd456a6"

    input:
    val(kraken_out)

    output:
    path("multi-krona.html")

    script:
    """
    echo "=========== Krona Plot ============"
    echo "Generating Krona plot from Kraken2 reports..."
    echo "Kraken files: ${kraken_out}"

    # Installing/updating krona database
    ktUpdateTaxonomy.sh

    # Running Krona on Kraken reports
    # -t is set to 7 as kraken reports are created with --report-minimizer-data
    ktImportTaxonomy -t 7 -m 3 -o multi-krona.html *_taxonomy.txt
    """
}
