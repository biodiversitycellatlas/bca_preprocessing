process SCIROCKET_SEQS_GATHER {
    tag "Gather demultiplexed sequencing for sample: ${meta.id}"

    input:
    tuple val(meta), file(scatter_dirs) from demux_scatter_ch

    output:
    tuple val(meta), file("demux_reads/${meta.id}") into gathered_seq_ch

    script:
    """
    mkdir -p demux_reads/${meta.id}
    
    # Combine QC pickles from all scatter chunks.
    python3 ${launchDir}/bin/scirocket_demux_gather.py \\
         --path_demux_scatter ${meta.id}/demux_reads_scatter/ \\
         --path_out demux_reads/${meta.id}/${meta.id}_qc.pickle

    # Combine discarded R1 and R2 reads as well as logs.
    find ./${meta.id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${meta.id}_R1_discarded.fastq.gz -print0 | xargs -0 cat > demux_reads/${meta.id}/${meta.id}_R1_discarded.fastq.gz
    find ./${meta.id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${meta.id}_R2_discarded.fastq.gz -print0 | xargs -0 cat > demux_reads/${meta.id}/${meta.id}_R2_discarded.fastq.gz
    find ./${meta.id}/demux_reads_scatter/ -maxdepth 2 -type f -name log_${meta.id}_discarded_reads.tsv.gz -print0 | xargs -0 cat > demux_reads/${meta.id}/log_${meta.id}_discarded_reads.tsv.gz

    # Combine and deduplicate whitelist files.
    find ./${meta.id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${meta.id}_whitelist_p7.txt -print0 | xargs -0 cat | sort -u > demux_reads/${meta.id}/${meta.id}_whitelist_p7.txt
    find ./${meta.id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${meta.id}_whitelist_p5.txt -print0 | xargs -0 cat | sort -u > demux_reads/${meta.id}/${meta.id}_whitelist_p5.txt
    find ./${meta.id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${meta.id}_whitelist_ligation.txt -print0 | xargs -0 cat | sort -u > demux_reads/${meta.id}/${meta.id}_whitelist_ligation.txt
    find ./${meta.id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${meta.id}_whitelist_rt.txt -print0 | xargs -0 cat | sort -u > demux_reads/${meta.id}/${meta.id}_whitelist_rt.txt
    """
}