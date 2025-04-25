process CONV_3CB_INDEX {
    publishDir "${params.output_dir}/demultiplex/demux_umitools/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("indexed_${fastqs[0]}")
    
    script:
    """
    echo "\n\n==================  CONV_3CB_INDEX =================="
    echo "Sample ID: ${meta}"
    echo "FASTQ files: ${fastqs}"

    python ${launchDir}/bin/conv_3cb_index.py \\
        --input_cDNA ${fastqs[0]} \\
        --output_cDNA indexed_${fastqs[0]} \\
        --workdir . \\
        --bcdir ${launchDir}/seq_techniques/bd_rhapsody
    """
}
