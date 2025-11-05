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
        // Create output channel for fastqs
        ch_demuxed_fastqs = Channel.create()

        // Demultiplex the fastq files based on the sample wells
        if (params.perform_demultiplexing && params.splitpipe_demultiplex_script == null) {
            log.info "Running Parse Biosciences demultiplexing using default script."
            PARSEBIO_CUSTOM_DEMUX(ch_samplesheet)
            ch_demuxed_fastqs = PARSEBIO_CUSTOM_DEMUX.out.splitted_files

        } else if (params.perform_demultiplexing && params.splitpipe_demultiplex_script != null) {
            log.info "Running Parse Biosciences demultiplexing using script: ${params.splitpipe_demultiplex_script}"
            PARSEBIO_PIPELINE_DEMUX(ch_samplesheet)
            ch_demuxed_fastqs = PARSEBIO_PIPELINE_DEMUX.out.splitted_files

        } else {
            log.info "Skipping Parse Biosciences demultiplexing as 'perform_demultiplexing' is set to false."
            ch_demuxed_fastqs = ch_samplesheet
        }

        // Only run Parse pipeline if the path is defined and exists
        if (params.splitpipe_installation && file(params.splitpipe_installation).exists()) {
            PARSEBIO_PIPELINE_MKREF()
            PARSEBIO_PIPELINE(ch_demuxed_fastqs, PARSEBIO_PIPELINE_MKREF.out)
        }

    emit:
        ch_demuxed_fastqs
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
