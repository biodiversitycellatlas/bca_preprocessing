process SAVE_RUN_CONFIG {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    label 'process_single'

    output:
    path("run_config_${params.trace_report_suffix}.txt")
    path("samplesheet_${params.trace_report_suffix}.csv")

    script:
    """
    # Copy the samplesheet 
    cp ${projectDir}/${params.input} "samplesheet_${params.trace_report_suffix}.csv"

    # Save the environment variables to a file
    cat <<EOF > run_config_${params.trace_report_suffix}.txt
    ${ params.keySet().sort()
        .collect { k -> "$k = ${params[k]}" }
        .join('\n') }
    EOF
    """
}
