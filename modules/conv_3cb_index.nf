process CONV_3CB_INDEX {
    publishDir "${params.output_dir}/demultiplex/demux_umitools/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val("${sample_id}"), path("indexed_*_R*.fastq.gz")
    
    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    """
    echo "\n\n==================  CONV_3CB_INDEX =================="
    echo "Sample ID: ${sample_id}"
    echo "FASTQ files: ${fastq_files}"

    python ${params.code_dir}/bin/conv_3cb_index.py \\
        --input_cDNA ${r2_fastq} \\
        --output_cDNA indexed_${sample_id}_R2.fastq.gz \\
        --workdir . \\
        --bcdir ${params.code_dir}/seq_techniques/bd_rhapsody
    """
}
