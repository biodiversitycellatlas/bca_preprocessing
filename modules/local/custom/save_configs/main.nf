process SAVE_RUN_CONFIG {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    label 'process_single2'
    tag "${input_file}"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    path("run_config_${params.trace_report_suffix}.txt"), emit: run_config
    path("samplesheet_${params.trace_report_suffix}.csv"), emit: samplesheet
    path "versions.yml", emit: versions

    script:
    """
    # Copy the samplesheet
    cp ${input_file} "samplesheet_${params.trace_report_suffix}.csv"

    # Save the environment variables to a file
    cat <<EOF > run_config_${params.trace_report_suffix}.txt
    ${ params.keySet().sort()
        .collect { k -> "$k = ${params[k]}" }
        .join('\n') }
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/^GNU bash, version //; s/ .*\$//')
END_VERSIONS
    """
}
