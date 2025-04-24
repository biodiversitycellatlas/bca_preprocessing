process SCIROCKET_SPL_GATHER {
    tag "Gather demultiplexed samples for sample: ${meta.id}"

    input:
    tuple val(meta), file(dummy) from gathered_seq_ch

    output:
    tuple val(meta),
          file("demux_reads/${meta.id}/${meta.id}_R1.fastq.gz"),
          file("demux_reads/${meta.id}/${meta.id}_R2.fastq.gz")
          into final_samples_ch

    script:
    """
    # Combine sample-specific FASTQ files from all demultiplexed splits.
    find ./${meta.id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${meta.id}_R1.fastq.gz -print0 | xargs -0 cat > demux_reads/${meta.id}/${meta.id}_R1.fastq.gz
    find ./${meta.id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${meta.id}_R2.fastq.gz -print0 | xargs -0 cat > demux_reads/${meta.id}/${meta.id}_R2.fastq.gz
    """
}