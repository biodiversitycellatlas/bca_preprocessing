process SCIROCKET_DEMUX {
    tag "Demultiplex sample: ${meta.id} chunk: ${chunk}"

    input:
    tuple val(meta), val(chunk), file(r1_chunk), file(r2_chunk) from demux_input_ch

    output:
    tuple val(meta), file("demux_reads_scatter/${meta.id}/${chunk}") into demux_scatter_ch

    script:
    """
    mkdir -p demux_reads_scatter/${meta.id}/${chunk}

    python3 ${launchDir}/bin/scirocket_demux_rocket.py \\
         --experiment_name ${meta.id} \\
         --samples ${params.path_samples} \\
         --barcodes ${params.path_barcodes} \\
         --r1 ${r1_chunk} --r2 ${r2_chunk} \\
         --out demux_reads_scatter/${meta.id}/${chunk} \\
         &> demultiplex_fastq_split_${meta.id}_${chunk}.log
    """
}
