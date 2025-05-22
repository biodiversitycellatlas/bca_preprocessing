process CONV_3CB_INDEX {
    publishDir "${params.output_dir}/demultiplex/demux_umitools/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)

    output:
    tuple val(meta), path("indexed_${fastq_cDNA}"), path(fastq_BC_UMI)
    
    script:
    // Get basename of the FASTQ files
    def fastq_cDNA_name = fastq_cDNA.toString().replaceAll(/.fastq.gz$/, '')
    def fastq_BC_UMI_name = fastq_BC_UMI.toString().replaceAll(/.fastq.gz$/, '')

    """
    echo "\n\n==================  CONV_3CB_INDEX =================="
    echo "Sample ID: ${meta}"
    echo "FASTQ files: ${fastqs}"

    python ${launchDir}/bin/conv_3cb_index.py \\
        --input_cDNA ${fastq_cDNA} \\
        --output_cDNA indexed_${fastq_cDNA_name} \\
        --workdir . \\
        --bcdir ${launchDir}/seq_techniques/bd_rhapsody
    """
}
