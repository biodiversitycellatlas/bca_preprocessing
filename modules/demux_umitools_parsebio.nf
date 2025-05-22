process DEMUX_UMITOOLS_PARSEBIO {
    publishDir "${params.output_dir}/demultiplex/demux_umitools/${meta.id}", mode: 'copy'
    tag "${meta.id}_${group}"
    debug true
    
    input:
    tuple val(meta), path(fastqs), val(group), val(wells)

    output:
    tuple val("${meta.id}_${group}"), path("demux_${fastq_cDNA_name}_${group}.fastq.gz"), path("demux_${fastq_BC_UMI_name}_${group}.fastq.gz"), emit: splitted_files

    script:
    // Get basename of the FASTQ files
    def fastq_cDNA_name = fastq_cDNA.toString().replaceAll(/.fastq.gz$/, '')
    def fastq_BC_UMI_name = fastq_BC_UMI.toString().replaceAll(/.fastq.gz$/, '')

    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${meta}"
    echo "First barcode path: ${params.barcode_umitools}"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "Wells: ${wells}"


    start_well=\$(echo ${wells} | cut -d'-' -f1)
    end_well=\$(echo ${wells} | cut -d'-' -f2)

    # Extract the row letter and start/end numbers
    row_letter="\${start_well:0:1}"
    start_num="\${start_well:1}"
    end_num="\${end_well:1}"

    rm -f ${group}_round1_whitelist.txt

    # Build the whitelist by looping over each well
    for ((i=start_num; i<=end_num; i++))
    do
        well="\${row_letter}\${i}"
        echo "> [\$well]"
        awk -F, -v well="\$well" '\$4 == well {print \$2}' "${params.barcode_umitools}" >> ${group}_round1_whitelist.txt
    done

    echo "whitelist: ${group}_round1_whitelist.txt "
    cat ${group}_round1_whitelist.txt

    umi_tools extract \\
        --extract-method=regex \\
        --bc-pattern2="^(?P<umi_1>.{10}){s<=1}.+(?P<cell_1>.{8}){s<=1}(?P<dummy>.)?\$" \\
        --stdin=${fastq_cDNA} \\
        --stdout=demux_${fastq_cDNA_name}_${group}_R1.fastq.gz \\
        --read2-in=${fastq_BC_UMI} \\
        --read2-out=demux_${fastq_BC_UMI_name}_${group}_R2.fastq.gz \\
        --whitelist=${group}_round1_whitelist.txt

    """
}

