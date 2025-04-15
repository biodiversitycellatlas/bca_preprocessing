process BDRHAP_PIPELINE {
    publishDir "${params.output_dir}/BDrhapsody_pipeline/${sample_id}", mode: 'copy', overwrite: false
    tag "${sample_id}_BDrhapsody"

    input:
    val sample_id
    path bd_ref_path

    script:
    """
    echo "\n\n===============  BD Rhapsody pipeline  ==============="
    echo "Sample: ${sample_id}"
    echo "BD Rhapsody reference files: ${bd_ref_path}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    # Add cwlref-runner to PATH
    export PATH=$PATH:\$CONDA_DEFAULT_ENV/env/lib/python3.13/site-packages

    cd ${params.bdrhap_pipeline_dir}

    # Define name of reference used during analysis
    basename_ref=\$(basename ${params.ref_gtf} .gtf)
    echo "basename: \${basename_ref}"

    cwl-runner \\
        --outdir ${params.output_dir}/BD_Rhapsody_pipeline/\${basename_ref} \\
        --singularity \\
        --reference-archive ${params.output_dir}/genome/BDrhap_reference/BD_Rhapsody_Reference_Files.tar.gz \\
        rhapsody_pipeline_2.2.1.cwl \\
        pipeline_inputs_${sample_id}_2.2.1.yml
    """
}