process CONV_3CB_INDEX {
    publishDir "${params.output_dir}/demultiplex/demux_umitools/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(fastq_files)

    output:
    tuple val(meta), path("indexed_*_R*.fastq.gz")
    
    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    """
    echo "\n\n==================  CONV_3CB_INDEX =================="
    echo "Sample ID: ${meta}"
    echo "FASTQ files: ${fastq_files}"

    python ${launchDir}/bin/conv_3cb_index.py \\
        --input_cDNA ${r2_fastq} \\
        --output_cDNA indexed_${meta.id}_R2.fastq.gz \\
        --workdir . \\
        --bcdir ${launchDir}/seq_techniques/bd_rhapsody
    """
}
