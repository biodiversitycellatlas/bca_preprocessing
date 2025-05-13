process SCIROCKET_DEMUX {
    tag "${meta.id}, ${fastqs}"
    debug true

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("demux_reads/${meta.id}_R{1,2}.fastq.gz"), emit: samples
    path("demux_reads/${meta.id}_whitelist_*"), emit: bc_whitelist

    script:
    // Retrieve barcode whitelist from conf/seqtech_parameters.config
    def bc_whitelist = params.seqtech_parameters[params.protocol].bc_whitelist

    """
    echo "\n\n==================  DEMULTIPLEXING FASTQ FILES  =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Sample ID: ${meta.id}"
    echo "Processing files: ${fastqs}"
    echo "Barcode whitelist: ${bc_whitelist}"
    echo "Samples file: ${launchDir}/${params.input}" 

    mkdir -p demux_reads/

    python3 ${launchDir}/bin/scirocket_demux_rocket.py \\
         --experiment_name ${meta.id} \\
         --samples ${launchDir}/${params.input} \\
         --barcodes ${bc_whitelist} \\
         --r1 ${fastqs[0]} --r2 ${fastqs[1]} \\
         --out demux_reads/ 
    """
}
