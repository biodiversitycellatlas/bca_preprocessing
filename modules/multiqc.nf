process MULTIQC {
    publishDir "${params.resDir}/multiqc", mode: 'copy'
    tag "all"
    
    input:
    val(trigger)

    output:
    file("multiqc_report.html")

    script:
    def rtitle = workflow.runName ? "--title \"${workflow.runName}\"" : ''
    def rfilename = workflow.runName ? "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    
    """    
    multiqc ${params.resDir} \\
        $rtitle \\
        $rfilename \\
        --config ${params.baseDir}/bin/multiqc_config.yaml \\
        --cl-config "report_header_info: '- Sequencing technology: ${params.seqTech}'"
    """
}
