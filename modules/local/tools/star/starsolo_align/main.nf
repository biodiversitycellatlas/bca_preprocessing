
// Retrieve one setting from custom parameters if set, otherwise from the protocol's block in conf/seqtech_parameters.config
def setting(Map seqtech, String key, Object fallback = null) {
    return params[key] != null ? params[key] : (seqtech[key] != null ? seqtech[key] : fallback)
}

process STARSOLO_ALIGN {
    publishDir "${params.outdir}/mapping_STARsolo/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high2'

    // Memory tracks the size of the FASTQs being read, overrides process_high2's flat assignments.
    // Coefficients live in params.dynamic_memory; remove the entry to fall back to the plain label.
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.STARSOLO_ALIGN,
        [fastq_cDNA, fastq_BC_UMI], task.attempt, 128, genome_index_files) }

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/htslib_samtools_seqspec_star_pruned:cef769c7e3b03dd0"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    val bc_whitelist
    val bam_only
    path genome_index_files

    output:
    tuple val(meta), path("*"),                                 emit: starsolo_files
    tuple val(meta), path("*_Log.final.out"),                   emit: log_final_file
    tuple val(meta), path("*_Log.out"),                         emit: log_file
    tuple val(meta), path("*_Solo.out"),                        emit: star_solodir
    tuple val(meta), path("*_Solo.out/GeneFull_Ex50pAS/raw"),   emit: genefull50_raw_dir, optional: true
    tuple val(meta), path("*_Solo.out/GeneFull_Ex50pAS/filtered"),      emit: genefull50_filtered_dir, optional: true
    tuple val(meta), path("*_Solo.out/GeneFull_Ex50pAS/Summary.csv"),   emit: summary_csv, optional: true
    tuple val(meta), path("*_Solo.out/GeneFull_Ex50pAS/CellReads.stats"), emit: cellreads_stats, optional: true
    tuple val(meta), path("*_Solo.out/GeneFull_Ex50pAS/UMIperCellSorted.txt"), emit: umi_per_cell, optional: true
    tuple val(meta), path("*_Solo.out/Velocyto/raw"),           emit: velocyto_raw_dir, optional: true
    tuple val(meta), path("*_Aligned.sortedByCoord.out.bam"),   emit: bam_file, optional: true
    path "versions.yml",                                        emit: versions

    script:
    // Retrieve settings from custom parameters if set, otherwise from conf/seqtech_parameters.config.
    def seqtech = params.seqtech_parameters[params.protocol]

    def star_soloTypestring    = setting(seqtech, 'star_soloTypestring')
    def star_soloCBmatchWLtype = setting(seqtech, 'star_soloCBmatchWLtype')
    def star_soloUMIfiltering  = setting(seqtech, 'star_soloUMIfiltering')
    def star_soloMultiMappers  = setting(seqtech, 'star_soloMultiMappers')
    def star_soloUMIdedup      = setting(seqtech, 'star_soloUMIdedup')
    def star_soloFeatures      = setting(seqtech, 'star_soloFeatures')
    def star_clipAdapterType   = setting(seqtech, 'star_clipAdapterType')
    def star_outFilterScoreMin = setting(seqtech, 'star_outFilterScoreMin')
    def star_outSAMunmapped    = setting(seqtech, 'star_outSAMunmapped')
    def star_outSAMattributes  = setting(seqtech, 'star_outSAMattributes')
    def star_solocellfilter    = setting(seqtech, 'star_solocellfilter')
    def star_extraargs         = setting(seqtech, 'star_extraargs', '')

    // Minimal STARsolo run with only one feature enabled for 'geneext_only' mode.
    def solo_features = bam_only ? ['Gene'] : star_soloFeatures.tokenize(' ').findAll { it }

    // params.perform_velocity adds the spliced/unspliced/ambiguous matrices used for RNA
    // velocity. STAR computes Velocyto from the Gene feature and exits without it, so Gene is
    // added back for a custom star_soloFeatures that left it out.
    if (params.perform_velocity && !bam_only) {
        if (!solo_features.contains('Gene')) { solo_features = ['Gene'] + solo_features }
        if (!solo_features.contains('Velocyto')) { solo_features << 'Velocyto' }
    }
    def star_soloFeatures_effective = solo_features.join(' ')

    def star_solocellfilter_effective  = bam_only ? 'None' : star_solocellfilter
    def star_soloCellReadStats         = bam_only ? 'None' : 'Standard'

    // Disable unmapped reads for 'geneext_only' mode.
    def star_outSAMunmapped_effective  = bam_only ? 'None' : star_outSAMunmapped

    // Convert empty bc_whitelist to None
    def safe_bc_whitelist = (bc_whitelist && bc_whitelist != "") ? bc_whitelist : 'None'

    // Detect CRAM or BAM by inspecting the first file in the list.
    def first_file     = fastq_cDNA instanceof List ? fastq_cDNA[0] : fastq_cDNA
    def is_bam_or_cram = first_file.name.endsWith(".bam") || first_file.name.endsWith(".cram")
    def read_files_command   = is_bam_or_cram ? "samtools view" : "pigz -dc -p ${task.cpus}"
    def input_files          = is_bam_or_cram
        ? (fastq_cDNA instanceof List ? fastq_cDNA.join(',') : "${fastq_cDNA}")
        : "${fastq_cDNA} ${fastq_BC_UMI}"

    // The barcode, UMI and gene tags are dropped when nothing will read them: either because
    // no BAM is written at all, or because the BAM exists only for GeneExt, which reads
    // coverage rather than tags. GX/GN/sM are in this list for the protocols that request them.
    def removableTags = ['CR', 'UR', 'CB', 'UB', 'GX', 'GN', 'sM']
    def star_outSAMattributes_effective = star_outSAMattributes
    if( params.star_generateBAM == false || bam_only ) {
        star_outSAMattributes_effective = star_outSAMattributes
            .split(/\s+/)
            .findAll { !(it in removableTags) }
            .join(' ')
    }

    // If star_generateBAM is false, do not output BAM (omit --outSAMtype), 'bam_only' forces it back on.
    def outSAMtype_option = (params.star_generateBAM || bam_only) ? '--outSAMtype BAM SortedByCoordinate' : '--outSAMtype None'

    // If cellfilter_method is "second_derivative" or "manual_cutoff", the filtered matrices are dropped and the raw matrices are used for cell calling.
    def drop_star_filtered = params.cellfilter_method in ["second_derivative", "manual_cutoff"]

    // Calculated in BcaResources.bamSortRam(), which returns a note and the byte count.
    // params.star_limitBAMsortRAM overrides whenever it is set to something non-zero.
    def bamsort = BcaResources.bamSortRam(params.star_limitBAMsortRAM, task.memory, genome_index_files)

    """
    echo "\n\n==============  MAPPING STARSOLO  ================"
    echo "Mapping sample ${meta.id} with STARsolo"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "Genome index directory: ${genome_index_files}"
    echo "Barcode whitelist: ${safe_bc_whitelist}"
    echo "Expected cells: ${meta.expected_cells}"
    echo "limitBAMsortRAM: ${bamsort.note}"
    echo "star_soloTypestring: ${star_soloTypestring}"
    echo "star_generateBAM: ${params.star_generateBAM}"
    echo "bam_only (GeneExt): ${bam_only}"
    echo "star_solocellfilter (effective): ${star_solocellfilter_effective}"
    echo "star_soloFeatures (effective): ${star_soloFeatures_effective}"
    echo "star_soloCellReadStats (effective): ${star_soloCellReadStats}"
    echo "star_outSAMunmapped (effective): ${star_outSAMunmapped_effective}"
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

    # Adjust soloCellFilter arguments. "None" takes no parameters of its own.
    SOLO_CELL_FILTER_ARGS=""
    if [[ "${star_solocellfilter_effective}" == "None" ]]; then
        SOLO_CELL_FILTER_ARGS="None"
    elif [[ "${star_solocellfilter_effective}" == "EmptyDrops_CR" ]]; then
        SOLO_CELL_FILTER_ARGS="${star_solocellfilter_effective} ${meta.expected_cells} 0.99 10 45000 90000 500 0.01 20000 0.01 10000"
    else
        SOLO_CELL_FILTER_ARGS="${star_solocellfilter_effective} ${meta.expected_cells} 0.99 10"
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
        --soloFeatures ${star_soloFeatures_effective} \\
        --soloCBwhitelist ${safe_bc_whitelist} \\
        --soloCellReadStats ${star_soloCellReadStats} \\
        --soloCellFilter \${SOLO_CELL_FILTER_ARGS} \\
        --clipAdapterType ${star_clipAdapterType} \\
        --outFilterScoreMin ${star_outFilterScoreMin} \\
        --outSAMunmapped ${star_outSAMunmapped_effective} \\
        --outSAMattributes ${star_outSAMattributes_effective} \\
        ${outSAMtype_option} \\
        --outFileNamePrefix ${meta.id}_ \\
        --genomeChrSetMitochondrial ${params.mt_contig} \\
        --limitBAMsortRAM ${bamsort.bytes} \\
        --soloStrand ${params.star_soloStrand} \\
        ${star_extraargs}

    if [[ "${drop_star_filtered}" == "true" ]]; then
        echo "cellfilter_method=${params.cellfilter_method}: removing STARsolo's own filtered matrices"
        rm -rf ${meta.id}_Solo.out/GeneFull_Ex50pAS/filtered
        rm -rf ${meta.id}_Solo.out/Gene/filtered
    fi

    # STAR insists on writing a raw matrix for the one feature it counted; nothing reads it
    # in this mode. Summary.csv and Barcodes.stats are kept as free QC on the alignment.
    if [[ "${bam_only}" == "true" ]]; then
        echo "bam_only: removing the unread Gene raw matrix, keeping Summary.csv and Barcodes.stats"
        rm -rf ${meta.id}_Solo.out/Gene/raw
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSIONS
    """
}
