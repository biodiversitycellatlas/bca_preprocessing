
/* Mapping using STAR (version 2.7.11b).             
 * The functions specify read 1 and read 2 from      
 * the channel, and otherwise state if not provided. 
 * The settings for STAR are set in the variable     
 * config_file, and is specific per sequencing tech. 
 * includes unmapped reads in the output.
*/

process MAPPING_STARSOLO { 
    publishDir "${params.output_dir}/mapping_STARsolo/${meta.id}", mode: 'copy', overwrite: false
    tag "${meta.id}_STARsolo"

    input:
    tuple val(meta), path(fastqs)
    path genome_index_files
    
    output:
    tuple val(meta), path("*")

    script:
    // Set default variables
    def bd_mem_arg = task.ext.args ?: ''

    // Retrieve barcode whitelist from conf/seqtech_parameters.config
    def bc_whitelist = params.seqtech_parameters[params.protocol].bc_whitelist

    // Retrieve starsolo settings from conf/seqtech_parameters.config
    def starsolo_settings = params.seqtech_parameters[params.protocol].starsolo
    def starsolo_args = starsolo_settings instanceof List ? starsolo_settings.join(' ') : starsolo_settings

    """
    echo "\n\n==============  MAPPING STARSOLO  ================"
    echo "Mapping sample ${meta.id} with STARsolo"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
    echo "Genome index directory: ${genome_index_files}"

    # In case the protocol does not exist and the user has not provided a seqspec file
    SOLO_ARGS=\"${starsolo_args}\"
    if [[ -n \"${params.seqspec_file}\" && \"${params.protocol}\" == *\"seqspec\"* ]]; 
    then
        # Retrieve input FASTQ files from seqspec file in starsolo format
        input_fastqs=\$(seqspec file -m rna -s read -f paired -k filename spec.yaml | awk '{print "-1 "\$1 " -2 "\$2}')

        # Barcode structure information from seqspec file
        bc_struct=\$(seqspec index -m rna -t starsolo -s file spec.yaml)
        SOLO_ARGS=\"\${bc_struct} \${SOLO_ARGS}\"
    else
        input_fastqs="-1 ${fastqs[0]} -2 ${fastqs[1]}"
    fi

    # Mapping step and generating count matrix using STAR
    STAR \\
        --runThreadN 8 \\
        --readFilesIn \${input_fastqs} \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${meta.id}_ \\
        --soloCellFilter CellRanger2.2 ${meta.expected_cells} 0.99 10 \\
        --soloCBwhitelist ${bc_whitelist} \\
        ${bd_mem_arg} \\
        \${SOLO_ARGS}
    """
}
