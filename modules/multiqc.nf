process MULTIQC {
    publishDir "${params.resDir}/multiqc", mode: 'copy'
    tag "all"
    
    input:
    val(trigger)

    output:
    file("multiqc_report.html")

    script:
    """    
    multiqc ${params.resDir} \\
        --config ${params.baseDir}/bin/multiqc_config.yml
    """
}
