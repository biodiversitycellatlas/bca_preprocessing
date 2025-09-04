process PARSEBIO_PIPELINE_DEMUX {
    publishDir "${params.outdir}/demultiplex/demux_spipe", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    
    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(input_file)

    output:
    tuple val(meta), path("*group*_R1*"), path("*group*_R2*"), path(input_file), emit: splitted_files

    script:
    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${meta}"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "Group: ${meta.id}"
    echo "Wells: ${meta.p5}" 

    # Run Parse Biosciences demultiplexing script
    fastq_sep_groups_v0.5.py \\
        --chemistry v3 \\
        --fq1 ${fastq_cDNA} \\
        --fq2 ${fastq_BC_UMI} \\
        --opath . \\
        --group ${meta.id} ${meta.p5} 
    """
}

