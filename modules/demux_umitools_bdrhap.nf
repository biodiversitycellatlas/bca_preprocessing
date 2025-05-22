process DEMUX_UMITOOLS_BDRHAP {
    publishDir "${params.output_dir}/demultiplex/demux_umitools/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true
    
    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)

    output:
    tuple val(meta), path("demux_${fastq_cDNA}"), path("demux_${fastq_BC_UMI}"), emit: splitted_files

    script:
    // Get basename of the FASTQ files
    def fastq_cDNA_name = fastq_cDNA.toString().replaceAll(/.fastq.gz$/, '')
    def fastq_BC_UMI_name = fastq_BC_UMI.toString().replaceAll(/.fastq.gz$/, '')

    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${meta}"
    echo "First barcode path: ${params.barcode_umitools}"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"

    umi_tools extract \\
        --extract-method=regex \\
        --bc-pattern="(?P<cell_1>.{9})GTGA(?P<cell_2>.{9})GACA(?P<cell_3>.{9})(?P<umi_1>.{8})" \\
        --stdin=${fastq_cDNA} \\
        --stdout=demux_${fastq_cDNA_name}.fastq.gz \\
        --read2-in=${fastq_BC_UMI} \\
        --read2-out=demux_${fastq_BC_UMI_name}.fastq.gz \\
        --log umi_extract.log

    """
}
