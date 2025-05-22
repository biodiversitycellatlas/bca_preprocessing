process PARSEBIO_CUSTOM_DEMUX {
    publishDir "${params.output_dir}/demultiplex/demux_custom/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true
    
    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)

    output:
    tuple val(meta), path("*group*_R1*"), path("*group*_R2*"), emit: splitted_files

    script:
    """
    echo "\n\n==================  Parse Biosciences: Custom Demultiplexing  =================="
    echo "Processing sample: ${meta}"
    echo "Fastq files: ${fastqs}"

    # Run Parse Biosciences demultiplexing script
    python ${launchDir}/bin/parsebio_custom_demux.py \\
        --sample_id ${meta.id} \\
        --fq1 ${fastq_cDNA} \\
        --fq2 ${fastq_BC_UMI} \\
        --whitelist ${launchDir}/seq_techniques/parse_biosciences/bc_data_n26_R1_v3_4.csv \\
        --group ${meta.group} ${meta.well} \\
        --output . \\
        --barcode_start 50 \\
        --barcode_end 58 
    """
}

