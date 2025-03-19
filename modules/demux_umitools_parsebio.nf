// ==================  DEMULTIPLEXING  ================== \\ 
// To split the fastq files of each library into separate \\
// fastq files for each fixation method, a script for     \\
// demultiplexing the reads is called. From a txt file    \\
// 'sample_wells', wells associated with each fixation    \\
// method are given. Based on the provided barcode file,  \\
// the samples are linked to the barcodes, and splitted   \\
// into n seperate fastq's. This step is repeated for     \\
// all libraries. 

process DEMUX_UMITOOLS_PARSEBIO {
    publishDir "${params.resDir}/demultiplex/demux_umitools/${sample_id}", mode: 'copy'
    tag "${sample_id}_${group}"
    debug true
    
    input:
    tuple val(sample_id), path(fastq_files), val(group), val(wells)

    output:
    tuple val("${sample_id}_${group}"), path("*_R*.fastq.gz"), emit: splitted_files

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${sample_id}"
    echo "First barcode path: ${params.barcode_umitools}"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
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
        --stdin=${r1_fastq} \\
        --stdout=demux_${sample_id}_${group}_R1.fastq.gz \\
        --read2-in=${r2_fastq} \\
        --read2-out=demux_${sample_id}_${group}_R2.fastq.gz \\
        --whitelist=${group}_round1_whitelist.txt

    """
}

