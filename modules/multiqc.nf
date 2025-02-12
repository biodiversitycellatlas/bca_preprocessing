// ==================  MULTI-QC  ================== \\ 
// A way of combining seperate FASTQC reports into  \\
// a single analysis. Provides an overview of the   \\
// data, including checks which of the files passed \\
// the quality metrics.                             \\

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
    echo "\n\n==================  Multi qc  =================="
    echo "Running MULTIQC"
    
    multiqc ${params.resDir} \\
        $rtitle \\
        $rfilename \\
        --config ${params.baseDir}/scripts/multiqc_config.yaml
        --cl-config "report_header_info: 'Sequencing technology: ${params.seqTech}'"
    """
}
