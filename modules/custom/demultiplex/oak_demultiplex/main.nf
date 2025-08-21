process OAK_DEMUX {
    publishDir "${params.output_dir}", mode: 'copy'
    tag "${meta.id}"
    label 'process_single'
    debug true
    
    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)

    output:
    tuple val(meta), path("demux_reads"), path("demux_reads/*.fastq.gz")

    when:

    script:
    """
    echo "\n\n==================  Demultiplex OAK data  =================="
    echo "Processing sample: ${meta}"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"

    echo "Demultiplexing reads..."

    # If i7 is provided, use it for demultiplexing first
    if [[ -n "${meta.p7}" ]]; then
        echo "Using i7 index: ${meta.p7}"
        bash ${moduleDir}/bin/oakseq_demux_i5_i7.sh \\
            ${meta.p7} \\
            i7 \\
            ${fastq_cDNA} \\
            ${fastq_BC_UMI} \\
            1
        echo "Demultiplexing i7 index completed."

        cDNA_file="demux_reads/cDNA_demux_i7.fastq.gz"
        BC_UMI_file="demux_reads/BC_UMI_demux_i7.fastq.gz"

    else
        echo "No i7 index provided, skipping i7 demultiplexing."
        cDNA_file="${fastq_cDNA}"
        BC_UMI_file="${fastq_BC_UMI}"
    fi

    # Always demultiplex i5
    echo "Using i5 index: ${meta.p5}"
    bash ${moduleDir}/bin/oakseq_demux_i5.sh \\
        ${meta.p5} \\
        i5 \\
        \${cDNA_file} \\
        \${BC_UMI_file} \\
        1

    echo "Demultiplexing i5 index completed."
    echo "Demultiplexing completed. Output files are in 'demux_reads' directory."
    """
}
