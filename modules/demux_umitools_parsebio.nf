process DEMUX_UMITOOLS_PARSEBIO {
    publishDir "${params.output_dir}/demultiplex/demux_umitools/${meta.id}", mode: 'copy'
    tag "${meta.id}_${group}"
    debug true
    
    input:
    tuple val(meta), path(fastqs), val(group), val(wells)

    output:
    tuple val("${meta.id}_${group}"), path("*_R*.fastq.gz"), emit: splitted_files

    script:
    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${meta}"
    echo "First barcode path: ${params.barcode_umitools}"
    echo "Fastqs: ${fastqs ?: 'Not provided'}"
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
        --stdin=${fastqs[0]} \\
        --stdout=demux_${meta.id}_${group}_R1.fastq.gz \\
        --read2-in=${fastqs[1]} \\
        --read2-out=demux_${meta.id}_${group}_R2.fastq.gz \\
        --whitelist=${group}_round1_whitelist.txt

    """
}

