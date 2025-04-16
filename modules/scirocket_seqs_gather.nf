process SCIROCKET_SEQS_GATHER {
    tag "Gather demultiplexed sequencing for sample: ${sample_id}"

    input:
    tuple val(sample_id), file(scatter_dirs) from demux_scatter_ch

    output:
    tuple val(sample_id), file("demux_reads/${sample_id}") into gathered_seq_ch

    script:
    """
    mkdir -p demux_reads/${sample_id}
    
    # Combine QC pickles from all scatter chunks.
    python3 ${params.code_dir}/bin/scirocket_demux_gather.py \\
         --path_demux_scatter ${sample_id}/demux_reads_scatter/ \\
         --path_out demux_reads/${sample_id}/${sample_id}_qc.pickle

    # Combine discarded R1 and R2 reads as well as logs.
    find ./${sample_id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${sample_id}_R1_discarded.fastq.gz -print0 | xargs -0 cat > demux_reads/${sample_id}/${sample_id}_R1_discarded.fastq.gz
    find ./${sample_id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${sample_id}_R2_discarded.fastq.gz -print0 | xargs -0 cat > demux_reads/${sample_id}/${sample_id}_R2_discarded.fastq.gz
    find ./${sample_id}/demux_reads_scatter/ -maxdepth 2 -type f -name log_${sample_id}_discarded_reads.tsv.gz -print0 | xargs -0 cat > demux_reads/${sample_id}/log_${sample_id}_discarded_reads.tsv.gz

    # Combine and deduplicate whitelist files.
    find ./${sample_id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${sample_id}_whitelist_p7.txt -print0 | xargs -0 cat | sort -u > demux_reads/${sample_id}/${sample_id}_whitelist_p7.txt
    find ./${sample_id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${sample_id}_whitelist_p5.txt -print0 | xargs -0 cat | sort -u > demux_reads/${sample_id}/${sample_id}_whitelist_p5.txt
    find ./${sample_id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${sample_id}_whitelist_ligation.txt -print0 | xargs -0 cat | sort -u > demux_reads/${sample_id}/${sample_id}_whitelist_ligation.txt
    find ./${sample_id}/demux_reads_scatter/ -maxdepth 2 -type f -name ${sample_id}_whitelist_rt.txt -print0 | xargs -0 cat | sort -u > demux_reads/${sample_id}/${sample_id}_whitelist_rt.txt
    """
}