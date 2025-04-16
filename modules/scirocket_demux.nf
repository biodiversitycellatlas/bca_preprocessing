process SCIROCKET_DEMUX {
    tag "Demultiplex sample: ${sample_id} chunk: ${chunk}"

    input:
    tuple val(sample_id), val(chunk), file(r1_chunk), file(r2_chunk) from demux_input_ch

    output:
    tuple val(sample_id), file("demux_reads_scatter/${sample_id}/${chunk}") into demux_scatter_ch

    script:
    """
    mkdir -p demux_reads_scatter/${sample_id}/${chunk}

    python3 ${params.code_dir}/bin/scirocket_demux_rocket.py \\
         --experiment_name ${sample_id} \\
         --samples ${params.path_samples} \\
         --barcodes ${params.path_barcodes} \\
         --r1 ${r1_chunk} --r2 ${r2_chunk} \\
         --out demux_reads_scatter/${sample_id}/${chunk} \\
         &> demultiplex_fastq_split_${sample_id}_${chunk}.log
    """
}
