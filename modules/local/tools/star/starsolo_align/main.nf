
process STARSOLO_ALIGN {
    publishDir "${params.outdir}/mapping_STARsolo/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/htslib_samtools_seqspec_star_pruned:cef769c7e3b03dd0"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    val bc_whitelist
    path genome_index_files

    output:
    tuple val(meta), path("*"),                                 emit: starsolo_files
    tuple val(meta), path("*_Log.final.out"),                   emit: log_final_file
    tuple val(meta), path("*_Log.out"),                         emit: log_file
    tuple val(meta), path("*_Solo.out"),                        emit: star_solodir
    tuple val(meta), path("*_Solo.out/GeneFull_Ex50pAS/raw"),   emit: genefull50_raw_dir
    tuple val(meta), path("*_Solo.out/GeneFull_Ex50pAS/filtered"),      emit: genefull50_filtered_dir
    tuple val(meta), path("*_Solo.out/GeneFull_Ex50pAS/Summary.csv"),   emit: summary_csv
    tuple val(meta), path("*_Solo.out/GeneFull_Ex50pAS/CellReads.stats"), emit: cellreads_stats
    tuple val(meta), path("*_Solo.out/GeneFull_Ex50pAS/UMIperCellSorted.txt"), emit: umi_per_cell
    tuple val(meta), path("*_Aligned.sortedByCoord.out.bam"),   emit: bam_file, optional: true

    script:
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

    // Convert empty bc_whitelist to None
    def safe_bc_whitelist = (bc_whitelist && bc_whitelist != "") ? bc_whitelist : 'None'

    // Detect CRAM or BAM by inspecting the first file in the list
    def first_file     = fastq_cDNA instanceof List ? fastq_cDNA[0] : fastq_cDNA
    def is_bam_or_cram = first_file.name.endsWith(".bam") || first_file.name.endsWith(".cram")
    def read_files_command   = is_bam_or_cram ? "samtools view" : "pigz -dc -p ${task.cpus}"
    def input_files          = is_bam_or_cram
        ? (fastq_cDNA instanceof List ? fastq_cDNA.join(',') : "${fastq_cDNA}")
        : "${fastq_cDNA} ${fastq_BC_UMI}"

    // If star_generateBAM is false, remove CR/UR/CB/UB tags from outSAMattributes
    def removableTags = ['CR', 'UR', 'CB', 'UB']
    def star_outSAMattributes_effective = star_outSAMattributes
    if( params.star_generateBAM == false) {
        star_outSAMattributes_effective = star_outSAMattributes
            .split(/\s+/)
            .findAll { !(it in removableTags) }
            .join(' ')
    }

    // If star_generateBAM is false, do not output BAM (omit --outSAMtype)
    def outSAMtype_option = params.star_generateBAM ? '--outSAMtype BAM SortedByCoordinate' : '--outSAMtype None'

    """
    echo "\n\n==============  MAPPING STARSOLO  ================"
    echo "Mapping sample ${meta.id} with STARsolo"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "Genome index directory: ${genome_index_files}"
    echo "Barcode whitelist: ${safe_bc_whitelist}"
    echo "Expected cells: ${meta.expected_cells}"
    echo "star_limitBAMsortRAM: ${params.star_limitBAMsortRAM}"
    echo "star_solocellfilter: ${star_solocellfilter}"
    echo "star_soloTypestring: ${star_soloTypestring}"
    echo "star_generateBAM: ${params.star_generateBAM}"
    echo "star_outSAMattributes (effective): ${star_outSAMattributes_effective}"
    echo "outSAMtype_option: ${outSAMtype_option}"

    if [[ -n \"${params.seqspec_file}\" && \"${params.protocol}\" == *\"seqspec\"* ]];
    then
        # Barcode structure information from seqspec file
        bc_struct=\$(seqspec index -m rna -t starsolo -s file ${params.seqspec_file})
        SOLO_TYPE_STRING="--soloType CB_UMI_Simple \${bc_struct}"
    else
        # Use predefined barcode structure information based on protocol
        SOLO_TYPE_STRING="${star_soloTypestring}"
    fi
    echo "SOLO_TYPE_STRING: \${SOLO_TYPE_STRING}"

    # Adjust soloCellFilter arguments
    SOLO_CELL_FILTER_ARGS=""
    if [[ "${star_solocellfilter}" == "EmptyDrops_CR" ]]; then
        SOLO_CELL_FILTER_ARGS="${star_solocellfilter} ${meta.expected_cells} 0.99 10 45000 90000 500 0.01 20000 0.01 10000"
    else
        SOLO_CELL_FILTER_ARGS="${star_solocellfilter} ${meta.expected_cells} 0.99 10"
    fi

    echo "SOLO_CELL_FILTER_ARGS: \${SOLO_CELL_FILTER_ARGS}"

    # Mapping step and generating count matrix using STAR
    STAR \\
        --runThreadN ${task.cpus} \\
        \${SOLO_TYPE_STRING} \\
        --readFilesIn ${input_files} \\
        --genomeDir ${genome_index_files} \\
        --readFilesCommand ${read_files_command} \\
        --soloCBmatchWLtype ${star_soloCBmatchWLtype} \\
        --soloUMIfiltering ${star_soloUMIfiltering} \\
        --soloMultiMappers ${star_soloMultiMappers} \\
        --soloUMIdedup ${star_soloUMIdedup} \\
        --soloFeatures ${star_soloFeatures} \\
        --soloCBwhitelist ${safe_bc_whitelist} \\
        --soloCellReadStats Standard \\
        --soloCellFilter \${SOLO_CELL_FILTER_ARGS} \\
        --clipAdapterType ${star_clipAdapterType} \\
        --outFilterScoreMin ${star_outFilterScoreMin} \\
        --outSAMunmapped ${star_outSAMunmapped} \\
        --outSAMattributes ${star_outSAMattributes_effective} \\
        ${outSAMtype_option} \\
        --outFileNamePrefix ${meta.id}_ \\
        --genomeChrSetMitochondrial ${params.mt_contig} \\
        --limitBAMsortRAM ${params.star_limitBAMsortRAM} \\
        --soloStrand ${params.star_soloStrand} \\
        ${star_extraargs}
    """
}
