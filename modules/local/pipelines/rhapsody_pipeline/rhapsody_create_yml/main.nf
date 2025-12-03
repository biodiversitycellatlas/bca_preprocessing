process BDRHAP_PIPELINE_YAML {
    label 'process_medium'

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path bd_ref_path

    output:
    tuple val(meta.id),                     emit: run_name
    path(bd_ref_path),                      emit: bd_ref_path
    path("rhapsody_input_${meta.id}.yml"),  emit: yaml_file

    script:
    """
    echo "\n\n===============  BD Rhapsody pipeline - create yaml  ==============="
    echo "Sample: ${meta.id}"
    echo "Fastq cDNA: ${fastq_cDNA}"
    echo "Fastq BC UMI: ${fastq_BC_UMI}"

    rhapsody_create_yaml.py \\
        --outprefix rhapsody_input_${meta.id} \\
        --fastq_cDNA ${fastq_cDNA} \\
        --fastq_BC_UMI ${fastq_BC_UMI} \\
        --star_ref ${bd_ref_path}
    """
}
