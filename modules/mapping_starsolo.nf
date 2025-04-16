
/* Mapping using STAR (version 2.7.11b).             
 * The functions specify read 1 and read 2 from      
 * the channel, and otherwise state if not provided. 
 * The settings for STAR are set in the variable     
 * config_file, and is specific per sequencing tech. 
 * includes unmapped reads in the output.
*/

process MAPPING_STARSOLO { 
    publishDir "${params.output_dir}/mapping_STARsolo/${sample_id}", mode: 'copy', overwrite: false
    tag "${sample_id}_STARsolo"

    input:
    tuple val(sample_id), path(fastq_files)
    path genome_index_files
    
    output:
    tuple val(sample_id), path("*")

    script:
    // Set default variables
    def bd_mem_arg = task.ext.args ?: ''

    def barcodeDir = "${params.code_dir}/seq_techniques/${params.protocol}/barcodes_R1.txt \
                        ${params.code_dir}/seq_techniques/${params.protocol}/barcodes_R2.txt \
                        ${params.code_dir}/seq_techniques/${params.protocol}/barcodes_R3.txt"
    def barcode_option = "--soloCBwhitelist ${barcodeDir.replaceAll('\n', '')}"

    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }
    def cDNA_read, CBUMI_read

    // Execute seqspec commands and capture output in Groovy variables
    if (params.protocol.toLowerCase().contains("bd_rhapsody")) {
        cDNA_read = r2_fastq
        CBUMI_read = r1_fastq
    
    } else if (params.protocol.toLowerCase().contains("oak_seq")) {
        cDNA_read = r2_fastq
        CBUMI_read = r1_fastq
        barcode_option = "--soloCBwhitelist ${params.code_dir}/seq_techniques/${params.protocol}/barcodes_R2.txt"
    
    } else {
        cDNA_read = r1_fastq
        CBUMI_read = r2_fastq
    }

    // If necesary, set the number of expected cells and cell filter for STARsolo
    if (!params.n_expected_cells || params.n_expected_cells == "") {
        n_expected_cells = 3000 // default STARsolo
    } else {
        n_expected_cells = params.n_expected_cells
    }
    def cell_filter = "--soloCellFilter CellRanger2.2 ${n_expected_cells} 0.99 10"

    """
    echo "\n\n==============  MAPPING STARSOLO  ================"
    echo "Mapping sample ${sample_id} with STARsolo"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
    echo "Genome index directory: ${genome_index_files}"
    echo "Barcodes: ${barcode_option}"

    # Read configuration file
    star_config="${params.code_dir}/seq_techniques/${params.protocol}/config_${params.protocol}_starsolo.txt"
    config_file=\$(cat \${star_config})

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
        ${cell_filter} \\
        \${config_file}
    """
}
