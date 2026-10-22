process PARSEBIO_CUSTOM_DEMUX {
    publishDir path: { "${params.outdir}/demultiplex/demux_custom/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930' :
        'community.wave.seqera.io/library/python:3.14.2--0562a0df3245213a' }"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path(whitelist, stageAs: 'splitwells_whitelist/*')

    output:
    tuple val(meta), path("*group*_R1*"), path("*group*_R2*"), path(fastq_indices), path(input_file), emit: splitted_files
    path "versions.yml", emit: versions

    script:
    """
    echo -e "\\n\\n==================  Parse Biosciences: Custom Demultiplexing  =================="
    echo "Processing sample: ${meta}"
    echo "Fastq files: ${fastq_cDNA}, ${fastq_BC_UMI}"
    echo "Group id: ${meta.id}, wells (rt): ${meta.rt}"

    # Run Parse Biosciences demultiplexing script
    parsebio_custom_demux.py \\
        --sample_id ${meta.id} \\
        --fq1 ${fastq_cDNA} \\
        --fq2 ${fastq_BC_UMI} \\
        --whitelist ${whitelist} \\
        --group ${meta.id} ${meta.rt} \\
        --output . \\
        --barcode_start 50 \\
        --barcode_end 58 \\
        --max_edit_dist 2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
