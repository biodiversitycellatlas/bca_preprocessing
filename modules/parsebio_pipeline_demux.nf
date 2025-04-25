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
    publishDir "${params.output_dir}/demultiplex/demux_spipe/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("*group*"), emit: splitted_files

    script:
    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${meta}"
    echo "FQ 1: ${fastqs[0] ?: 'Not provided'}"
    echo "FQ 2: ${fastqs[1] ?: 'Not provided'}"

    # Run Parse Biosciences demultiplexing script
    python ${launchDir}/bin/fastq_sep_groups_v0.5.py \\
        --chemistry v3 \\
        --fq1 ${fastqs[0]} \\
        --fq2 ${fastqs[1]} \\
        --opath . \\
        --group ${meta.group} ${meta.well} 
    """
}

