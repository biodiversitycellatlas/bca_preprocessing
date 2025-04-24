process MULTIQC {
    publishDir "${params.output_dir}/multiqc", mode: 'copy'
    tag "all"
    
    input:
    val(trigger)

    output:
    file("multiqc_report.html")

    script:
    """    
    multiqc ${params.output_dir} \\
        --config ${launchDir}/bin/multiqc_config.yml
    """
}
