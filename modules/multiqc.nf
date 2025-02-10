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
    """
    echo "\n\n==================  Multi qc  =================="
    echo "Running MULTIQC"
    ls -l

    multiqc ${params.resDir} \\
        --dirs \\
        --dirs-depth 3 \\
        -v \\
        -f \\
        -m cutadapt -m star -m fastqc -m featureCounts
    """
}
