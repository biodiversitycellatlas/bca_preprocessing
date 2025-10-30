process PARSEBIO_CUSTOM_DEMUX {
    publishDir "${params.outdir}/demultiplex/demux_custom/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'


    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(input_file)

    output:
    tuple val(meta), path("*group*_R1*"), path("*group*_R2*"), path(input_file), emit: splitted_files

    script:
    """
    echo "\n\n==================  Parse Biosciences: Custom Demultiplexing  =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Processing sample: ${meta}"
    echo "Fastq files: ${fastq_cDNA}, ${fastq_BC_UMI}"

    # Run Parse Biosciences demultiplexing script
    parsebio_custom_demux.py \\
        --sample_id ${meta.id} \\
        --fq1 ${fastq_cDNA} \\
        --fq2 ${fastq_BC_UMI} \\
        --whitelist ${params.bc_whitelist_parse_splitwells} \\
        --group ${meta.id} ${meta.p5} \\
        --output . \\
        --barcode_start 50 \\
        --barcode_end 58
    """
}
