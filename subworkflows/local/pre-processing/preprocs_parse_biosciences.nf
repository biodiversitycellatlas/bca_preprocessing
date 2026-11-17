//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PARSEBIO_CUSTOM_DEMUX } from '../../../modules/local/custom/demultiplex/parsebio_demultiplex/main'
include { PARSEBIO_COLLAPSE_BC1 } from '../../../modules/local/custom/demultiplex/parsebio_collapse_bc1/main'
include { PARSEBIO_PIPELINE_DEMUX } from '../../../modules/local/pipelines/split-pipe/split-pipe_demux/main'
include { PARSEBIO_PIPELINE_MKREF } from '../../../modules/local/pipelines/split-pipe/split-pipe_mkref/main'
include { PARSEBIO_PIPELINE } from '../../../modules/local/pipelines/split-pipe/split-pipe_all/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR PARSE BIOSCIENCES DATA
        Based on the sample wells, the fastq sequences must be split
        based on the first barcode round. Performs either custom- and
        commercial pre-processing (using Parse Biosciences pipeline v1.3.1)
        enabling comparison between the methods and validation of steps.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow parse_workflow {
    take:
        ch_samplesheet

    main:
        def ch_versions = channel.empty()
        def ch_splitpipe_stats = channel.empty()

        // Demultiplex the fastq files based on the sample wells
        if (params.perform_demultiplexing && params.splitpipe_demultiplex_script == null) {
            log.info "Running Parse Biosciences demultiplexing using default script."
            if (!params.bc_whitelist_parse_splitwells) {
                error("The Parse Biosciences demultiplexing splits the reads on the wells in 'bc_whitelist_parse_splitwells', which is not set.")
            }
            PARSEBIO_CUSTOM_DEMUX(ch_samplesheet, file(params.bc_whitelist_parse_splitwells, checkIfExists: true))
            demux_samplesheet = PARSEBIO_CUSTOM_DEMUX.out.splitted_files
            ch_versions = ch_versions.mix(PARSEBIO_CUSTOM_DEMUX.out.versions)

        } else if (params.perform_demultiplexing && params.splitpipe_demultiplex_script != null) {
            log.info "Running Parse Biosciences demultiplexing using script: ${params.splitpipe_demultiplex_script}"
            PARSEBIO_PIPELINE_DEMUX(ch_samplesheet)
            demux_samplesheet = PARSEBIO_PIPELINE_DEMUX.out.splitted_files
            ch_versions = ch_versions.mix(PARSEBIO_PIPELINE_DEMUX.out.versions)

        } else {
            log.info "Skipping Parse Biosciences demultiplexing as 'perform_demultiplexing' is set to false."
            demux_samplesheet = ch_samplesheet
        }

        // Each round-1 well holds an oligo-dT (T), a random-hexamer (R) and in some kits
        // targeted (X) RT primers, each with its own bc1. split-pipe names cells by well, so
        // bc1 is rewritten to the well's T barcode before mapping: otherwise each cell is
        // split over its T and R barcodes. split-pipe itself keeps the original reads.
        def ch_mapping_input = demux_samplesheet
        def ch_bc1_whitelist = channel.empty()
        if (params.parse_collapse_bc1) {
            if (!params.bc_whitelist_parse_splitwells) {
                error("Merging the round-1 barcodes of each well ('parse_collapse_bc1') reads the wells from 'bc_whitelist_parse_splitwells', which is not set. Set it, or set 'parse_collapse_bc1 = false'.")
            }
            PARSEBIO_COLLAPSE_BC1(
                demux_samplesheet.map { meta, _fastq_cDNA, fastq_BC_UMI, _fastq_indices, _input_file -> [meta, fastq_BC_UMI] },
                file(params.bc_whitelist_parse_splitwells, checkIfExists: true)
            )
            ch_versions = ch_versions.mix(PARSEBIO_COLLAPSE_BC1.out.versions)

            // Swap the rewritten barcode read into each sample's tuple
            ch_mapping_input = demux_samplesheet
                .map { meta, fastq_cDNA, _fastq_BC_UMI, fastq_indices, input_file -> [meta.id, meta, fastq_cDNA, fastq_indices, input_file] }
                .join(PARSEBIO_COLLAPSE_BC1.out.reads.map { meta, collapsed -> [meta.id, collapsed] })
                .map { _id, meta, fastq_cDNA, fastq_indices, input_file, collapsed -> [meta, fastq_cDNA, collapsed, fastq_indices, input_file] }

            // The T barcode of every well: the round-1 whitelist for the mappers (the same for every sample)
            ch_bc1_whitelist = PARSEBIO_COLLAPSE_BC1.out.bc1_whitelist.first()
        }

        // Only run Parse pipeline if the path is defined and exists
        if (params.splitpipe_installation && file(params.splitpipe_installation).exists()) {
            PARSEBIO_PIPELINE_MKREF(file(params.ref_fasta, checkIfExists: true), file(params.ref_gtf, checkIfExists: true))

            // The reference process has single-value inputs only, so it runs once and its output is reused by every sample
            PARSEBIO_PIPELINE(demux_samplesheet, PARSEBIO_PIPELINE_MKREF.out.reference)
            ch_versions = ch_versions.mix(PARSEBIO_PIPELINE_MKREF.out.versions, PARSEBIO_PIPELINE.out.versions)

            // The two split-pipe reports the mapping statistics are read from
            ch_splitpipe_stats = PARSEBIO_PIPELINE.out.sample_stats.mix(PARSEBIO_PIPELINE.out.agg_summary)
        }

    emit:
        data_output     = ch_mapping_input
        bc1_whitelist   = ch_bc1_whitelist
        splitpipe_stats = ch_splitpipe_stats
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
