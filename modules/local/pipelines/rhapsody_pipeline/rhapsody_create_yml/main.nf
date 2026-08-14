process BDRHAP_PIPELINE_YAML {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path bd_ref_path

    output:
    tuple val(meta), path(bd_ref_path), path("rhapsody_input_${meta.id}.yml"), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file), emit: pipeline_input
    path "versions.yml",                    emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
