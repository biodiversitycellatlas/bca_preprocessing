process BDRHAP_PIPELINE {
    publishDir "${params.outdir}/BDrhapsody_pipeline/${meta.id}", mode: 'copy'
    tag "${meta.id}_BDrhapsody"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"

    input:
    val meta
    path bd_ref_path

    script:
    """
    echo "\n\n===============  BD Rhapsody pipeline  ==============="
    echo "Sample: ${meta}"
    echo "BD Rhapsody reference files: ${bd_ref_path}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    # Add cwlref-runner to PATH
    export PATH=$PATH:\$CONDA_DEFAULT_ENV/env/lib/python3.13/site-packages

    cd ${params.rhapsody_installation}

    # Define name of reference used during analysis
    basename_ref=\$(basename ${params.ref_gtf} .gtf)
    echo "basename: \${basename_ref}"

    cwl-runner \\
        --outdir ${params.outdir}/BD_Rhapsody_pipeline/\${basename_ref} \\
        --singularity \\
        --reference-archive ${params.outdir}/genome/BDrhap_reference/BD_Rhapsody_Reference_Files.tar.gz \\
        rhapsody_pipeline_2.2.1.cwl \\
        pipeline_inputs_${meta.id}_2.2.1.yml
    """
}
