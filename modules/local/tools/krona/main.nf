process KRONA {
    publishDir "${params.outdir}/summary_results", mode: 'copy'
    label 'process_single'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.8--pl5262hdfd78af_2' :
        'biocontainers/krona:2.8--pl5262hdfd78af_2' }"

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
