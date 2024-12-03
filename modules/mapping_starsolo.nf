// ==============  MAPPING STARSOLO  =============== \\ 
// Mapping using STAR (version 2.7.11b).             \\
// The functions specify read 1 and read 2 from      \\
// the channel, and otherwise state if not provided. \\
// The settings for STAR are set in the variable     \\
// config_file, and is specific per sequencing tech. \\

process MAPPING_STARSOLO {    
    publishDir "${params.resDir}/mapping_STARsolo_${config_name}/${sample_id}", mode: 'symlink'
    debug true
    label 'big_mem'
    tag "${fastq_files}_${config_name}"

    input:
    tuple val(sample_id), path(fastq_files)
    path genome_index_files
    path star_config
    val config_name
    
    output:
    tuple val(sample_id), val(config_name), path("*")

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }
    """
    echo "\n\n==============  MAPPING STARSOLO  ================"
    echo "Mapping sample ${sample_id} with STARsolo"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
    echo "Genome index directory: ${genome_index_files}"

    # Read configuration file
    config_file=\$(cat ${star_config})

    # Mapping step and generating count matrix using STAR
    STAR \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${r1_fastq.baseName}_ \\
        --readFilesIn ${r1_fastq} ${r2_fastq} \\
        --soloCBwhitelist ${params.barcodeDir} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes CR UR CB UB \\
        \${config_file} 
    """
}
