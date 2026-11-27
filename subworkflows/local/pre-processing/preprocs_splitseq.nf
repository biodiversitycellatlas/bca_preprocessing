//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SPLITSEQ_WELLS        } from '../../../modules/local/custom/demultiplex/splitseq_wells/main'
include { SPLITSEQ_DEMUX        } from '../../../modules/local/custom/demultiplex/splitseq_demultiplex/main'
include { SPLITSEQ_COLLAPSE_BC1 } from '../../../modules/local/custom/demultiplex/splitseq_collapse_bc1/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR SPLiT-seq DATA (Rosenberg et al. 2018)
        Read 2 holds UMI, round-3, round-2 and round-1 (RT) barcodes; each RT well
        holds an oligo-dT and a random-hexamer primer with barcodes of their own.
        - Samples with 'rt' wells are split out of their sub-library's reads
        - The round-1 barcodes of each well are merged into its oligo-dT barcode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow splitseq_workflow {
    take:
        ch_samplesheet  // [meta, fastq_cDNA, fastq_BC_UMI, fastq_indices, input_file]
        ch_whitelist    // the 'bc_whitelist' files: one list for all rounds, or bc1, bc2, bc3

    main:
        def ch_versions = channel.empty()
        def seqtech   = params.seqtech_parameters[params.protocol]
        def bc1_start = seqtech.splitseq_bc1_start
        def bc1_end   = seqtech.splitseq_bc1_end

        // Plain barcode lists for the mappers and the round-1 well CSV
        def wells_given = params.splitseq_wells ? file(params.splitseq_wells, checkIfExists: true) : []
        SPLITSEQ_WELLS(ch_whitelist, wells_given)
        ch_versions = ch_versions.mix(SPLITSEQ_WELLS.out.versions)
        def ch_wells = SPLITSEQ_WELLS.out.wells

        // A sample with 'rt' wells is split out of its sub-library's reads; one without
        // takes every read (one sub-library, one sample)
        def ch_split = ch_samplesheet.branch { row ->
            demux: params.perform_demultiplexing && (row[0].rt?.toString()?.trim() as boolean)
            whole: true
        }
        if (!params.perform_demultiplexing) {
            log.info "SPLiT-seq: 'perform_demultiplexing' is false, so every sample takes all reads of its FASTQs, whatever its 'rt' wells."
        }
        SPLITSEQ_DEMUX(ch_split.demux, ch_wells, bc1_start, bc1_end)
        ch_versions = ch_versions.mix(SPLITSEQ_DEMUX.out.versions)
        def ch_reads = ch_split.whole.mix(SPLITSEQ_DEMUX.out.splitted_files)

        // Each RT well holds an oligo-dT and a random-hexamer barcode: rewrite bc1 to the
        // well's oligo-dT barcode before mapping, or every cell is split in two
        def ch_mapping_input = ch_reads
        def ch_bc1_whitelist = SPLITSEQ_WELLS.out.bc1
        if (params.splitseq_collapse_bc1) {
            SPLITSEQ_COLLAPSE_BC1(
                ch_reads.map { meta, _fastq_cDNA, fastq_BC_UMI, _fastq_indices, _input_file -> [meta, fastq_BC_UMI] },
                ch_wells, bc1_start, bc1_end
            )
            ch_versions = ch_versions.mix(SPLITSEQ_COLLAPSE_BC1.out.versions)

            // Swap the rewritten barcode read into each sample's tuple
            ch_mapping_input = ch_reads
                .map { meta, fastq_cDNA, _fastq_BC_UMI, fastq_indices, input_file -> [meta.id, meta, fastq_cDNA, fastq_indices, input_file] }
                .join(SPLITSEQ_COLLAPSE_BC1.out.reads.map { meta, collapsed -> [meta.id, collapsed] })
                .map { _id, meta, fastq_cDNA, fastq_indices, input_file, collapsed -> [meta, fastq_cDNA, collapsed, fastq_indices, input_file] }

            // The oligo-dT barcode of every well (the same for every sample)
            ch_bc1_whitelist = SPLITSEQ_COLLAPSE_BC1.out.bc1_whitelist.first()
        }

        // One list per segment, in the order of 'soloCBposition' and the alevin geometry
        def ch_bc_whitelist = ch_bc1_whitelist
            .combine(SPLITSEQ_WELLS.out.bc2)
            .combine(SPLITSEQ_WELLS.out.bc3)
            .map { bc1, bc2, bc3 -> [bc1, bc2, bc3] }
            .first()

    emit:
        data_output  = ch_mapping_input
        bc_whitelist = ch_bc_whitelist
        versions     = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
