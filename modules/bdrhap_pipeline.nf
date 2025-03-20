process BDRHAP_PIPELINE {
    publishDir "${params.resDir}/BD_Rhapsody_pipeline", mode: 'copy'

    input:
    file("BD_Rhapsody_Reference_Files.tar.gz")

    script:
    """
    echo "\n\n===============  BD Rhapsody pipeline  ==============="
    cd ${params.bdrhap_pipeline_dir}

    cwl-runner \\
        --outdir . \\
        --singularity \\
        rhapsody_pipeline_2.2.1.cwl \\
        pipeline_inputs_BCA005_2.2.1.yml
    """
}