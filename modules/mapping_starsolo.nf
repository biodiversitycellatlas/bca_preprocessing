
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

    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    // Retrieve barcode whitelist from conf/seqtech_parameters.config
    def bc_whitelist = params.seqtech_parameters[params.protocol].bc_whitelist

    // Retrieve starsolo settings from conf/seqtech_parameters.config
    def starsolo_settings = params.seqtech_parameters[params.protocol].starsolo
    def starsolo_args = starsolo_settings instanceof List ? starsolo_settings.join(' ') : starsolo_settings

    // Define which FASTQ file contains cDNA and which are the barcodes/UMIs
    def cDNA_read, CBUMI_read
    if (params.protocol.contains("bd_rhapsody") || params.protocol.contains("10x") || params.protocol.contains("oak_seq") || params.protocol.contains("sci_rna_seq3")) {
        cDNA_read = r2_fastq
        CBUMI_read = r1_fastq
    } else {
        cDNA_read = r1_fastq
        CBUMI_read = r2_fastq
    }

    """
    echo "\n\n==============  MAPPING STARSOLO  ================"
    echo "Mapping sample ${sample_id} with STARsolo"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
    echo "Genome index directory: ${genome_index_files}"

    # Mapping step and generating count matrix using STAR
    STAR \\
        --runThreadN 8 \\
        --readFilesIn ${cDNA_read} ${CBUMI_read} \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${sample_id}_ \\
        --soloCellFilter CellRanger2.2 ${params.n_expected_cells} 0.99 10 \\
        --soloCBwhitelist ${bc_whitelist} \\
        ${bd_mem_arg} \\
        ${starsolo_args}
    """
}
