process BDRHAP_PIPELINE {
    publishDir "${params.resDir}/BD_Rhapsody_pipeline/\${basename_ref}", mode: 'copy', overwrite: false

    input:
    path(BD_ref_path)

    script:
    """
    echo "\n\n===============  BD Rhapsody pipeline  ==============="
    echo ${BD_ref_path}

    cd ${params.bdrhap_pipeline_dir}

    # Define name of reference genome
    basename_ref=\$(basename ${params.ref_star_gtf})

    cwl-runner \\
        --outdir . \\
        --singularity \\
        rhapsody_pipeline_2.2.1.cwl \\
        pipeline_inputs_BCA005_2.2.1.yml
    """
}