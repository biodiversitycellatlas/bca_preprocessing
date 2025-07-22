
/* Mapping using STAR (version 2.7.11b).             
 * The functions specify read 1 and read 2 from      
 * the channel, and otherwise state if not provided. 
 * The settings for STAR are set in the variable     
 * config_file, and is specific per sequencing tech. 
 * includes unmapped reads in the output.
*/

process MAPPING_STARSOLO { 
    publishDir "${params.output_dir}/mapping_STARsolo/${meta.id}", mode: 'copy'
    tag "${meta.id}_STARsolo"
    debug true

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)
    val bc_whitelist
    path genome_index_files
    
    output:
    tuple val(meta), path("*")

    script:
    // Set default variables
    def bd_mem_arg = task.ext.args ?: ''

    // Retrieve starsolo settings from conf/seqtech_parameters.config
    def starsolo_settings = params.seqtech_parameters[params.protocol].starsolo
    def starsolo_args = starsolo_settings instanceof List ? starsolo_settings.join(' ') : starsolo_settings

    """
    echo "\n\n==============  MAPPING STARSOLO  ================"
    echo "Mapping sample ${meta.id} with STARsolo"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "Genome index directory: ${genome_index_files}"
    echo "Barcode whitelist: ${bc_whitelist}"
    echo "Expected cells: ${meta.expected_cells}"

    # In case the protocol does not exist and the user has not provided a seqspec file
    SOLO_ARGS=\"${starsolo_args}\"
    if [[ -n \"${params.seqspec_file}\" && \"${params.protocol}\" == *\"seqspec\"* ]]; 
    then
        # Barcode structure information from seqspec file
        bc_struct=\$(seqspec index -m rna -t starsolo -s file spec.yaml)
        SOLO_ARGS=\"\${bc_struct} \${SOLO_ARGS}\"
    fi

    echo "SOLO_ARGS: \${SOLO_ARGS}"
    echo "First 10 lines of FASTQ cDNA:"
    zcat ${fastq_cDNA} | head -n 10
    echo "First 10 lines of FASTQ BC & UMI:"
    zcat ${fastq_BC_UMI} | head -n 10

    # Mapping step and generating count matrix using STAR
    STAR \\
        --runThreadN 8 \\
        --readFilesIn ${fastq_cDNA} ${fastq_BC_UMI} \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand "pigz -dc -p ${task.cpus}" \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${meta.id}_ \\
        --soloCellFilter CellRanger2.2 ${meta.expected_cells} 0.99 10 \\
        --soloCBwhitelist ${bc_whitelist} \\
        ${bd_mem_arg} \\
        \${SOLO_ARGS}
    """
}
