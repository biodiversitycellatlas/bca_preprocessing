process BDRHAP_PIPELINE {
    publishDir "${params.outdir}/BDrhapsody_pipeline/", mode: 'copy'
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    val run_name
    path bd_ref_path
    path input_yaml

    script:
    """
    echo "\n\n===============  BD Rhapsody pipeline  ==============="
    echo "Run name: ${run_name}"
    echo "BD Rhapsody reference files: ${bd_ref_path}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    cwltool \\
        --outdir ${params.outdir}/BD_Rhapsody_pipeline/${run_name} \\
        --singularity \\
        ${params.rhapsody_installation}/rhapsody_pipeline_2.2.1.cwl \\
        ${input_yaml}
    """
}
