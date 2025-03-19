process DEMUX_UMITOOLS_BDRHAP {
    publishDir "${params.resDir}/demultiplex/demux_umitools/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true
    
    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val("${sample_id}"), path("*_R*.fastq.gz"), emit: splitted_files

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${sample_id}"
    echo "First barcode path: ${params.barcode_umitools}"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"

    umi_tools extract \\
        --extract-method=regex \\
        --bc-pattern="(?P<cell_1>.{9})GTGA(?P<cell_2>.{9})GACA(?P<cell_3>.{9})(?P<umi_1>.{8})" \\
        --stdin=${r1_fastq} \\
        --stdout=demux_${sample_id}_R1.fastq.gz \\
        --read2-in=${r2_fastq} \\
        --read2-out=demux_${sample_id}_R2.fastq.gz \\
        --log umi_extract.log

    """
}
