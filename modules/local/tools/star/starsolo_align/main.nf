
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

    // Retrieve settings from custom parameters if set, otherwise from conf/seqtech_parameters.config
    def star_soloTypestring = params.star_soloTypestring ?: params.seqtech_parameters[params.protocol].star_soloTypestring
    def star_soloCBmatchWLtype = params.star_soloCBmatchWLtype ?: params.seqtech_parameters[params.protocol].star_soloCBmatchWLtype
    def star_soloUMIfiltering = params.star_soloUMIfiltering ?: params.seqtech_parameters[params.protocol].star_soloUMIfiltering
    def star_soloMultiMappers = params.star_soloMultiMappers ?: params.seqtech_parameters[params.protocol].star_soloMultiMappers
    def star_soloUMIdedup = params.star_soloUMIdedup ?: params.seqtech_parameters[params.protocol].star_soloUMIdedup
    def star_soloFeatures = params.star_soloFeatures ?: params.seqtech_parameters[params.protocol].star_soloFeatures
    def star_clipAdapterType = params.star_clipAdapterType ?: params.seqtech_parameters[params.protocol].star_clipAdapterType
    def star_outFilterScoreMin = params.star_outFilterScoreMin ?: params.seqtech_parameters[params.protocol].star_outFilterScoreMin
    def star_outSAMunmapped = params.star_outSAMunmapped ?: params.seqtech_parameters[params.protocol].star_outSAMunmapped
    def star_outSAMattributes = params.star_outSAMattributes ?: params.seqtech_parameters[params.protocol].star_outSAMattributes
    def star_solocellfilter = params.star_solocellfilter ?: params.seqtech_parameters[params.protocol].star_solocellfilter
    def star_extraargs = params.star_extraargs ?: params.seqtech_parameters[params.protocol].star_extraargs

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
    if [[ "${star_solocellfilter}" == "EmptyDrops_CR" ]]; then
        SOLO_CELL_FILTER_ARGS="--soloCellFilter ${star_solocellfilter} ${meta.expected_cells} 0.99 10 45000 90000 500 0.01 20000 0.01 10000"
    else
        SOLO_CELL_FILTER_ARGS="--soloCellFilter ${star_solocellfilter} ${meta.expected_cells} 0.99 10"
    fi

    # Mapping step and generating count matrix using STAR
    STAR \\
        --runThreadN 8 \\
        --readFilesIn ${fastq_cDNA} ${fastq_BC_UMI} \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand "pigz -dc -p ${task.cpus}" \\
        --soloCBmatchWLtype ${star_soloCBmatchWLtype} \\
        --soloUMIfiltering ${star_soloUMIfiltering} \\
        --soloMultiMappers ${star_soloMultiMappers} \\
        --soloUMIdedup ${star_soloUMIdedup} \\
        --soloFeatures ${star_soloFeatures} \\
        --soloCBwhitelist ${bc_whitelist} \\
        --soloCellReadStats Standard \\
        --soloCellFilter ${SOLO_CELL_FILTER_ARGS} \\
        --clipAdapterType ${star_clipAdapterType} \\
        --outFilterScoreMin ${star_outFilterScoreMin} \\
        --outSAMunmapped ${star_outSAMunmapped} \\
        --outSAMattributes ${star_outSAMattributes} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${meta.id}_ \\
        --genomeChrSetMitochondrial ${params.mt_contig} \\
        ${limitBAMsortRAM} \\
        ${star_extraargs}
    """
}
