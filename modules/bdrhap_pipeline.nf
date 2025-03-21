process BDRHAP_PIPELINE {
    publishDir "${params.resDir}/BD_Rhapsody_pipeline/\${basename_ref}", mode: 'copy', overwrite: false

    input:
    val sample_id
    path bd_ref_path

    script:
    """
    echo "\n\n===============  BD Rhapsody pipeline  ==============="
    echo "Sample: ${sample_id}"
    echo "BD Rhapsody reference files: ${bd_ref_path}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    cd ${params.bdrhap_pipeline_dir}

    # Define name of reference used during analysis
    basename_ref=\$(basename ${params.ref_star_gtf} .gtf)
    echo "basename: \${basename_ref}"

    cwl-runner \\
        --outdir . \\
        --singularity \\
        rhapsody_pipeline_2.2.1.cwl \\
        pipeline_inputs_${sample_id}_2.2.1.yml
    """
}