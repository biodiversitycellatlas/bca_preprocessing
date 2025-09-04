
process STARSOLO_ALIGN { 
    publishDir "${params.outdir}/mapping_STARsolo/${meta.id}", mode: 'copy'
    tag "${meta.id}_STARsolo"
    label 'process_high'
    

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(input_file)
    val bc_whitelist
    path genome_index_files
    
    output:
    tuple val(meta), path("*")

    script:
    // Set default variables
    def limitBAMsortRAM = (params.protocol == 'bd_rhapsody' || params.protocol == '10xv3') ? '--limitBAMsortRAM 50000000000' : '' 

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
    echo "limitBAMsortRAM: ${limitBAMsortRAM}"

    # In case the protocol does not exist and the user has not provided a seqspec file
    SOLO_ARGS=\"${starsolo_args}\"
    if [[ -n \"${params.seqspec_file}\" && \"${params.protocol}\" == *\"seqspec\"* ]]; 
    then
        # Barcode structure information from seqspec file
        bc_struct=\$(seqspec index -m rna -t starsolo -s file spec.yaml)
        SOLO_ARGS=\"\${bc_struct} \${SOLO_ARGS}\"
    fi

    echo "SOLO_ARGS: \${SOLO_ARGS}"

    # Adjust soloCellFilter arguments
    SOLO_CELL_FILTER_ARGS=""
    if [[ "${params.star_solocellfilter}" == "EmptyDrops_CR" ]]; then
        SOLO_CELL_FILTER_ARGS="--soloCellFilter ${params.star_solocellfilter} ${meta.expected_cells} 0.99 10 45000 90000 500 0.01 20000 0.01 10000"
    else
        SOLO_CELL_FILTER_ARGS="--soloCellFilter ${params.star_solocellfilter} ${meta.expected_cells} 0.99 10"
    fi

    # Mapping step and generating count matrix using STAR
    STAR \\
        --runThreadN 8 \\
        --readFilesIn ${fastq_cDNA} ${fastq_BC_UMI} \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand "pigz -dc -p ${task.cpus}" \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${meta.id}_ \\
        \${SOLO_CELL_FILTER_ARGS} \\
        --soloCBwhitelist ${bc_whitelist} \\
        --soloCellReadStats Standard \\
        --genomeChrSetMitochondrial ${params.mt_contig} \\
        ${limitBAMsortRAM} \\
        \${SOLO_ARGS}
    """
}
