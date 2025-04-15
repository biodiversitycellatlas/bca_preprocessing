// ==================  DEMUX_SPIPE  ================== \\ 
// To split the fastq files of each library into separate \\
// fastq files for each fixation method, a script for     \\
// demultiplexing the reads is called. From a txt file    \\
// 'sample_wells', wells associated with each fixation    \\
// method are given. Based on the provided barcode file,  \\
// the samples are linked to the barcodes, and splitted   \\
// into n seperate fastq's. This step is repeated for     \\
// all libraries. 

process PARSEBIO_PIPELINE_DEMUX {
    publishDir "${params.output_dir}/demultiplex/demux_spipe/${sample_id}", mode: 'copy'
    tag "${sample_id}_${group}"
    
    input:
    tuple val(sample_id), path(fastq_files), val(group), val(wells)

    output:
    tuple val("${sample_id}_${group}"), path("*group*"), emit: splitted_files

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }
    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${sample_id}"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"

    # Run Parse Biosciences demultiplexing script
    python ${params.code_dir}/bin/fastq_sep_groups_v0.5.py \\
        --chemistry v3 \\
        --fq1 ${r1_fastq} \\
        --fq2 ${r2_fastq} \\
        --opath . \\
        --group ${group} ${wells} 
    """
}

