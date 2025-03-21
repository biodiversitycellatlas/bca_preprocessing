
/* Mapping using STAR (version 2.7.11b).             
 * The functions specify read 1 and read 2 from      
 * the channel, and otherwise state if not provided. 
 * The settings for STAR are set in the variable     
 * config_file, and is specific per sequencing tech. 
*/

process MAPPING_STARSOLO { 
    publishDir "${params.resDir}/mapping_STARsolo/\${basename_ref}_incl_unmapped/${sample_id}", mode: 'copy', overwrite: false
    tag "${sample_id}_STARsolo"

    input:
    tuple val(sample_id), path(fastq_files)
    path genome_index_files
    
    output:
    tuple val(sample_id), path("*")

    script:
    def bd_mem_arg   = task.ext.args ?: '' 
    def barcode_option = "--soloCBwhitelist ${params.barcodeDir.replaceAll('\n', '')}" 

    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    // If seqTech is "bd_rhapsody", then cDNA = R2 and CB/UMI = R1, else by default cDNA = R1 and CB/UMI = R2
    def cDNA_read
    def CBUMI_read
    if (params.seqTech.toLowerCase().contains("bd_rhapsody")) {
        cDNA_read = r2_fastq
        CBUMI_read = r1_fastq
    } else if (params.seqTech.toLowerCase().contains("oak_seq")) {
        cDNA_read = r2_fastq
        CBUMI_read = r1_fastq
        barcode_option = "--soloCBwhitelist ${params.baseDir}/seq_techniques/${params.seqTech}/barcodes_R2.txt"
    } else {
        cDNA_read = r1_fastq
        CBUMI_read = r2_fastq
    }

    """
    echo "\n\n==============  MAPPING STARSOLO ================"
    echo "Mapping sample ${sample_id} with STARsolo"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
    echo "Genome index directory: ${genome_index_files}"
    echo "Barcodes: ${barcode_option}"

    # Read configuration file
    config_file=\$(cat ${params.star_config})

    # Define name of reference genome
    basename_ref=\$(basename ${params.ref_star_gtf} .gtf)

    # Mapping step and generating count matrix using STAR
    STAR \\
        --runThreadN 8 \\
        --readFilesIn ${cDNA_read} ${CBUMI_read} \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand zcat \\
        ${barcode_option} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS nM CR UR CB UB \\
        --soloMultiMappers EM \\
        --outFileNamePrefix ${sample_id}_ \\
        ${bd_mem_arg} \\
        \${config_file} 
    """
}
