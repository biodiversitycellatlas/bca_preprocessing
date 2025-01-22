// ==================  MULTI-QC  ================== \\ 
// A way of combining seperate FASTQC reports into  \\
// a single analysis. Provides an overview of the   \\
// data, including checks which of the files passed \\
// the quality metrics.                             \\

process MULTIQC {
    publishDir "${params.resDir}/fastqc", mode: 'copy'
    tag "all"
    
    input:
    path('*')

    output:
    file("multiqc_report.html")

    script:
    """
    echo "\n\n==================  Multi qc  =================="
    echo "Running MULTIQC"
    multiqc .
    """
}
