process SCIROCKET_SPL_GATHER {
    tag "Gather demultiplexed samples for sample: ${meta.id}"
    debug true

    input:
    tuple val(meta), file(dummy)

    output:
    tuple val(meta), path("demux_reads/demux_R{1,2}.fastq.gz")

    script:
    """
    echo "\n\n==================  GATHER DEMULTIPLEXED SAMPLES  =================="
    echo "Sample ID: ${meta}"
    echo "Processing files: ${dummy}"
    echo "Gathering demultiplexed samples from ${dummy}"
    echo "Output directory: demux_reads/"

    # Combine sample-specific FASTQ files from all demultiplexed splits.
    find ./demux_reads_scatter/ -maxdepth 2 -type f -name *_R1.fastq.gz -print0 | xargs -0 cat > demux_reads/demux_R1.fastq.gz
    find ./demux_reads_scatter/ -maxdepth 2 -type f -name *_R2.fastq.gz -print0 | xargs -0 cat > demux_reads/demux_R2.fastq.gz
    """
}