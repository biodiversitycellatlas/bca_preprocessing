process SCIROCKET_SEQS_GATHER {
    debug true

    input:
    path(samples_discarded)
    path(bc_whitelists)

    output:
    path("seqs_gather/"),                   emit: all_output
    path("seqs_gather/whitelist_*"),        emit: bc_whitelist

    script:
    """
    echo "\n\n==================  GATHER DEMULTIPLEXED SEQUENCING DATA  =================="
    echo "Discarded sample files: ${samples_discarded}"
    echo "Barcode whitelist files: ${bc_whitelists}"
    echo "Output directory: seqs_gather/"

    mkdir -p seqs_gather/
    
    # Combine QC pickles from all scatter chunks.
    # python3 ${launchDir}/bin/scirocket_demux_gather.py \\
    #     --path_demux_scatter demux_reads/ \\
    #     --path_out seqs_gather/qc.pickle

    # Combine discarded R1 and R2 reads as well as logs.
    find . -type f -name *_R1_discarded.fastq.gz -print0 | xargs -0 cat > seqs_gather/R1_discarded.fastq.gz
    find . -type f -name *_R2_discarded.fastq.gz -print0 | xargs -0 cat > seqs_gather/R2_discarded.fastq.gz

    # Combine and deduplicate whitelist files.
    find . -type f -name *_whitelist_p7.txt -print0 | xargs -0 cat | sort -u > seqs_gather/whitelist_p7.txt
    find . -type f -name *_whitelist_p5.txt -print0 | xargs -0 cat | sort -u > seqs_gather/whitelist_p5.txt
    find . -type f -name *_whitelist_ligation.txt -print0 | xargs -0 cat | sort -u > seqs_gather/whitelist_ligation.txt
    find . -type f -name *_whitelist_rt.txt -print0 | xargs -0 cat | sort -u > seqs_gather/whitelist_rt.txt
    """
}