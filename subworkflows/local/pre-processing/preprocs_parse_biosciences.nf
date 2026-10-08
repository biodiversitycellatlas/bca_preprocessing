//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PARSEBIO_CUSTOM_DEMUX } from '../../../modules/local/custom/demultiplex/parsebio_demultiplex/main'
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
            PARSEBIO_CUSTOM_DEMUX(ch_samplesheet)
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

        // Only run Parse pipeline if the path is defined and exists
        if (params.splitpipe_installation && file(params.splitpipe_installation).exists()) {
            PARSEBIO_PIPELINE_MKREF()

            // Use .first() to reuse the reference output for all split samples
            PARSEBIO_PIPELINE(demux_samplesheet, PARSEBIO_PIPELINE_MKREF.out.reference.first())
            ch_versions = ch_versions.mix(PARSEBIO_PIPELINE_MKREF.out.versions, PARSEBIO_PIPELINE.out.versions)

            // The two split-pipe reports the mapping statistics are read from
            ch_splitpipe_stats = PARSEBIO_PIPELINE.out.sample_stats.mix(PARSEBIO_PIPELINE.out.agg_summary)
        }

    emit:
        data_output     = demux_samplesheet
        splitpipe_stats = ch_splitpipe_stats
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
