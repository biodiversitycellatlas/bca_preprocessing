process DEMUX_UMITOOLS_BDRHAP {
    publishDir "${params.output_dir}/demultiplex/demux_umitools/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true
    
    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("*_R*.fastq.gz"), emit: splitted_files

    script:
    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${meta}"
    echo "First barcode path: ${params.barcode_umitools}"
    echo "Fastq files: ${fastqs ?: 'Not provided'}"

    umi_tools extract \\
        --extract-method=regex \\
        --bc-pattern="(?P<cell_1>.{9})GTGA(?P<cell_2>.{9})GACA(?P<cell_3>.{9})(?P<umi_1>.{8})" \\
        --stdin=${fastqs[0]} \\
        --stdout=demux_${meta.id}_R1.fastq.gz \\
        --read2-in=${fastqs[1]} \\
        --read2-out=demux_${meta.id}_R2.fastq.gz \\
        --log umi_extract.log

    """
}
