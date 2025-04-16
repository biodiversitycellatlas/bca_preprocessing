process SCIROCKET_SPL_GATHER {
    tag "Gather demultiplexed samples for sample: ${sample_id}"

    input:
    tuple val(sample_id), file(dummy) from gathered_seq_ch

    output:
    tuple val(sample_id),
          file("demux_reads/${sample_id}/${sample_id}_R1.fastq.gz"),
          file("demux_reads/${sample_id}/${sample_id}_R2.fastq.gz")
          into final_samples_ch

    script:
    """
    # Combine sample-specific FASTQ files from all demultiplexed splits.
    find ./${sample_id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${sample_id}_R1.fastq.gz -print0 | xargs -0 cat > demux_reads/${sample_id}/${sample_id}_R1.fastq.gz
    find ./${sample_id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${sample_id}_R2.fastq.gz -print0 | xargs -0 cat > demux_reads/${sample_id}/${sample_id}_R2.fastq.gz
    """
}