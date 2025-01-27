// ==============  MAPPING STARSOLO  =============== \\ 
// Mapping using STAR (version 2.7.11b).             \\
// The functions specify read 1 and read 2 from      \\
// the channel, and otherwise state if not provided. \\
// The settings for STAR are set in the variable     \\
// config_file, and is specific per sequencing tech. \\

process MAPPING_STARSOLO { 
    publishDir "${params.resDir}/mapping_STARsolo/mapping_STARsolo_${config_name}/${sample_id}", mode: 'copy', overwrite: false
    tag "${sample_id}_${config_name}"

    input:
    tuple val(sample_id), path(fastq_files)
    path genome_index_files
    path star_config
    val barcode_whitelist_input
    val config_name
    
    output:
    tuple val(sample_id), val(config_name), path("*")

    script:
    def bd_mem_arg   = task.ext.args ?: ''          // If ext.args is defined assign it to bd_mem_arg

    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    def barcode_option = barcode_whitelist_input ? "--soloCBwhitelist ${barcode_whitelist_input.replaceAll('\n', '')}" : ""

    // If config_name contains "bd_rhapsody", then cDNA = R2 and CB/UMI = R1
    // Else by default cDNA = R1 and CB/UMI = R2
    def cDNA_read
    def CBUMI_read
    if (params.seqTech.toLowerCase().contains("bd_rhapsody")) {
        cDNA_read = r2_fastq
        CBUMI_read = r1_fastq
    } else {
        cDNA_read = r1_fastq
        CBUMI_read = r2_fastq
    }

    """
    echo "\n\n==============  MAPPING STARSOLO ${config_name} ================"
    echo "Mapping sample ${sample_id} with STARsolo"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
    echo "Genome index directory: ${genome_index_files}"
    echo "Barcodes: ${barcode_whitelist_input}"
    echo "Barcode input STAR: ${barcode_option}"

    # Read configuration file
    config_file=\$(cat ${star_config})

    # Mapping step and generating count matrix using STAR
    STAR \\
        --runThreadN 40 \\
        --readFilesIn ${cDNA_read} ${CBUMI_read} \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand zcat \\
        ${barcode_option} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS nM CR UR CB UB \\
        --soloMultiMappers EM \\
        ${bd_mem_arg} \\
        \${config_file} 

    """
}
