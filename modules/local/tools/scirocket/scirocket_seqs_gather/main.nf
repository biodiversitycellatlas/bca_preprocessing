process SCIROCKET_SEQS_GATHER {
    publishDir "${params.outdir}", mode: 'copy'
    label 'process_single'


    input:
    path samples_discarded_R1
    path samples_discarded_R2
    path bc_whitelists_p5
    path bc_whitelists_p7
    path bc_whitelists_ligation
    path bc_whitelists_rt

    output:
    tuple path('seqs_gather/whitelist_p7.txt'), path('seqs_gather/whitelist_p5.txt'), path('seqs_gather/whitelist_ligation.txt'), path('seqs_gather/whitelist_rt.txt'), emit: bc_whitelist

    script:
    """
    echo "\n\n==================  GATHER DEMULTIPLEXED SEQUENCING DATA  =================="
    echo "Output directory: seqs_gather/"

    mkdir -p seqs_gather/

    # Combine discarded R1 and R2 reads as well as logs.
    cat ${samples_discarded_R1.join(' ')} > seqs_gather/R1_discarded.fastq.gz
    cat ${samples_discarded_R2.join(' ')} > seqs_gather/R2_discarded.fastq.gz

    # Combine and deduplicate whitelist files.
    cat ${bc_whitelists_p5.join(' ')} | sort -u > seqs_gather/whitelist_p5.txt
    cat ${bc_whitelists_p7.join(' ')} | sort -u > seqs_gather/whitelist_p7.txt
    cat ${bc_whitelists_ligation.join(' ')} | sort -u > seqs_gather/whitelist_ligation.txt
    cat ${bc_whitelists_rt.join(' ')} | sort -u > seqs_gather/whitelist_rt.txt
    """
}
