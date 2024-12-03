// ==================  DEMULTIPLEXING  ================== \\ 
// To split the fastq files of each library into separate \\
// fastq files for each fixation method, a script for     \\
// demultiplexing the reads is called. From a txt file    \\
// 'sample_wells', wells associated with each fixation    \\
// method are given. Based on the provided barcode file,  \\
// the samples are linked to the barcodes, and splitted   \\
// into n seperate fastq's. This step is repeated for     \\
// all libraries. 

process DEMULTIPLEX {
    publishDir "${params.resDir}/demux_fastq", mode: 'symlink'
    tag "${sample_id}"
    label 'big_mem'
    debug true
    
    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: splitted_files

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }
    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${sample_id}"
    echo "First barcode path: ${params.barcodeDemux}"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"

    # Run the demux code
    bash ${params.codeDir}/scripts/split_parse_data.sh ${params.resDir} ${sample_id} ${r1_fastq} ${r2_fastq} ${params.barcodeDemux}
    """
}

