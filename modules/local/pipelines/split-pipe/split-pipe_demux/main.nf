process PARSEBIO_PIPELINE_DEMUX {
    publishDir path: { "${params.outdir}/demultiplex/demux_spipe" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/numpy_pandas_python:b78d8279027c3c38' :
        'community.wave.seqera.io/library/numpy_pandas_python:e63d9be2ea3d0055' }"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    tuple val(meta), path("*group*_R1*"), path("*group*_R2*"), path(fastq_indices), path(input_file), emit: splitted_files
    path "versions.yml", emit: versions

    script:
    // Retrieve kit specification for split-pipe from custom parameters if set, otherwise from conf/seqtech_parameters.config
    def kit_specification = params.splitpipe_kit ? "--kit ${params.splitpipe_kit} --kit_score_skip" : ''

    // The script belongs to the split-pipe installation, so it is run where it is rather than
    // staged; made absolute, since the task runs in its own work directory
    def demux_script = file(params.splitpipe_demultiplex_script).toString()

    """
    echo -e "\\n\\n==================  split-pipe Demultiplex  =================="
    echo "split-pipe script: ${demux_script}"
    echo "Processing sample: ${meta}"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "Group: ${meta.id}"
    echo "Wells: ${meta.rt}"

    # Run Parse Biosciences demultiplexing script
    python ${demux_script} \\
        --chemistry v3 \\
        --fq1 ${fastq_cDNA} \\
        --fq2 ${fastq_BC_UMI} \\
        --opath . \\
        --group ${meta.id} ${meta.rt} \\
        ${kit_specification}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
