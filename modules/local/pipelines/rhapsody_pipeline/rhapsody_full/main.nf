process BDRHAP_PIPELINE {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/BDrhapsody_pipeline/" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cwltool_python:d7b534cc8e4f1511' :
        'community.wave.seqera.io/library/cwltool_python:910578d904c1d231' }"

    input:
    tuple val(meta), path(bd_ref_path), path(input_yaml), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    path "versions.yml",                 emit: versions

    script:
    def run_name = meta.id
    // The workflow belongs to the BD Rhapsody installation and references files next to it, so it
    // is run where it is rather than staged; made absolute, since the task runs in its own work directory
    def rhapsody_cwl = "${file(params.rhapsody_installation)}/rhapsody_pipeline_2.2.1.cwl"
    """
    echo -e "\\n\\n===============  BD Rhapsody pipeline  ==============="
    echo "Run name: ${run_name}"
    echo "BD Rhapsody reference files: ${bd_ref_path}"

    # Written into the task directory and published from there, like every other module
    cwltool \\
        --outdir ${run_name} \\
        --singularity \\
        ${rhapsody_cwl} \\
        ${input_yaml}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cwltool: \$(cwltool --version 2>&1 | sed 's/^.*cwltool //')
    END_VERSIONS
    """
}
