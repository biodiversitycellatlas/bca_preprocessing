process SCIROCKET_SPL_GATHER {
    publishDir "${params.output_dir}/demux_reads/spl_gather", mode: 'copy'
    tag "all-samples"
    debug true

    input:
    path r1_list
    path r2_list

    output:
    tuple val("all"), path("spl_gather/all_R1.fastq.gz"), path("spl_gather/all_R2.fastq.gz")

    script:
    """
    echo "\n\n==================  GATHER DEMULTIPLEXED SAMPLES  =================="
    echo "Gathering ${r1_list.size()} R1 files and ${r2_list.size()} R2 files"
    echo "Output directory: spl_gather/"

    ls
    mkdir -p spl_gather/

    # Combine sample-specific FASTQ files from all demultiplexed splits.
    cat ${r1_list.join(' ')} > spl_gather/all_R1.fastq.gz
    cat ${r2_list.join(' ')} > spl_gather/all_R2.fastq.gz
    """
}