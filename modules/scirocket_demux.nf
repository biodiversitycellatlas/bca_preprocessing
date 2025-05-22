process SCIROCKET_DEMUX {
    publishDir "${params.output_dir}", mode: 'copy'
    tag "${meta.id}, ${fastq_cDNA}, ${fastq_BC_UMI}"
    debug true

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)

    output:
    path("demux_reads/${meta.id}_R{1,2}.fastq.gz"),       emit: samples
    path("demux_reads/${meta.id}_R{1,2}_discarded*"),     emit: samples_discarded
    path("demux_reads/${meta.id}_whitelist_*"),           emit: bc_whitelists

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
