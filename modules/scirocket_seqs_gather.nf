process SCIROCKET_SEQS_GATHER {
    tag "Gather demultiplexed sequencing for sample: ${meta.id}"
    debug true

    input:
    tuple val(meta), path(scatter_dirs)

    output:
    tuple val(meta), path("demux_reads/")

    script:
    """
    echo "\n\n==================  GATHER DEMULTIPLEXED SEQUENCING DATA  =================="
    echo "Sample ID: ${meta}"
    echo "Processing files: ${scatter_dirs}"
    echo "Gathering demultiplexed sequencing data from ${scatter_dirs}"
    echo "Output directory: demux_reads/"

    mkdir -p demux_reads
    
    # Combine QC pickles from all scatter chunks.
    python3 ${launchDir}/bin/scirocket_demux_gather.py \\
         --path_demux_scatter demux_reads_scatter/ \\
         --path_out demux_reads/qc.pickle

    # Combine discarded R1 and R2 reads as well as logs.
    find ./demux_reads_scatter/ -maxdepth 2 -type f -name R1_discarded.fastq.gz -print0 | xargs -0 cat > demux_reads/R1_discarded.fastq.gz
    find ./demux_reads_scatter/ -maxdepth 2 -type f -name R2_discarded.fastq.gz -print0 | xargs -0 cat > demux_reads/R2_discarded.fastq.gz
    find ./demux_reads_scatter/ -maxdepth 2 -type f -name log_discarded_reads.tsv.gz -print0 | xargs -0 cat > demux_reads/log_discarded_reads.tsv.gz

    # Combine and deduplicate whitelist files.
    find ./demux_reads_scatter/ -maxdepth 2 -type f -name whitelist_p7.txt -print0 | xargs -0 cat | sort -u > demux_reads/whitelist_p7.txt
    find ./demux_reads_scatter/ -maxdepth 2 -type f -name whitelist_p5.txt -print0 | xargs -0 cat | sort -u > demux_reads/whitelist_p5.txt
    find ./demux_reads_scatter/ -maxdepth 2 -type f -name whitelist_ligation.txt -print0 | xargs -0 cat | sort -u > demux_reads/whitelist_ligation.txt
    find ./demux_reads_scatter/ -maxdepth 2 -type f -name whitelist_rt.txt -print0 | xargs -0 cat | sort -u > demux_reads/whitelist_rt.txt
    """
}