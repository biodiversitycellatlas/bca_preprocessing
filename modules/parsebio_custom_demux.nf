process PARSEBIO_CUSTOM_DEMUX {
    publishDir "${params.resDir}/demultiplex/demux_custom/${sample_id}", mode: 'copy'
    tag "${sample_id}_${group}"
    debug true
    
    input:
    tuple val(sample_id), path(fastq_files), val(group), val(wells)

    output:
    tuple val("${sample_id}_${group}"), path("*group*"), emit: splitted_files

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }
    """
    echo "\n\n==================  Parse Biosciences: Custom Demultiplexing  =================="
    echo "Processing sample: ${sample_id}"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
    echo "Group: ${group}"
    echo "Wells: ${wells}"

    # Run Parse Biosciences demultiplexing script
    python ${params.baseDir}/bin/parsebio_custom_demux.py \\
        --sample_id ${sample_id} \\
        --fq1 ${r1_fastq} \\
        --fq2 ${r2_fastq} \\
        --whitelist ${params.baseDir}/seq_techniques/parse_biosciences/bc_data_n26_R1_v3_4.csv \\
        --group ${group} ${wells} \\
        --output . \\
        --barcode_start 50 \\
        --barcode_end 58 
    """
}

