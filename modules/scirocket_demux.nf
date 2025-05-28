process SCIROCKET_DEMUX {
    publishDir "${params.output_dir}", mode: 'copy'
    tag "${meta.id}, ${fastq_cDNA}, ${fastq_BC_UMI}"
    debug true

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)

    output:
    path("demux_reads/${meta.id}_R1.fastq.gz"),       emit: samples_R1
    path("demux_reads/${meta.id}_R2.fastq.gz"),       emit: samples_R2
    path("demux_reads/${meta.id}_R1_discarded*"),     emit: samples_discarded_R1
    path("demux_reads/${meta.id}_R2_discarded*"),     emit: samples_discarded_R2
    path("demux_reads/${meta.id}_whitelist_p5*"),           emit: bc_whitelists_p5
    path("demux_reads/${meta.id}_whitelist_p7*"),           emit: bc_whitelists_p7
    path("demux_reads/${meta.id}_whitelist_ligation*"),     emit: bc_whitelists_ligation
    path("demux_reads/${meta.id}_whitelist_rt*"),           emit: bc_whitelists_rt           

    script:
    // Retrieve barcode whitelist from conf/seqtech_parameters.config
    def bc_whitelist = params.seqtech_parameters[params.protocol].bc_whitelist

    // Get basename of the FASTQ files
    def fastq_cDNA_name = fastq_cDNA.toString().replaceAll(/.fastq.gz$/, '')
    def fastq_BC_UMI_name = fastq_BC_UMI.toString().replaceAll(/.fastq.gz$/, '')

    """
    echo "\n\n==================  DEMULTIPLEXING FASTQ FILES  =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Sample ID: ${meta.id}"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "Barcode whitelist: ${bc_whitelist}"
    echo "Samples file: ${launchDir}/${params.input}" 

    mkdir -p demux_reads/

    python3 ${launchDir}/bin/scirocket_demux_rocket.py \\
         --experiment_name ${meta.id} \\
         --samples ${launchDir}/${params.input} \\
         --barcodes ${bc_whitelist} \\
         --r1 ${fastq_BC_UMI} --r2 ${fastq_cDNA} \\
         --out demux_reads/ 
    """
}
