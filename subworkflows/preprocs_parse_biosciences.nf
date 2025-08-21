//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PARSEBIO_PIPELINE_DEMUX } from '../modules/pipelines/split-pipe/split-pipe_demux/main'
include { PARSEBIO_PIPELINE_MKREF } from '../modules/pipelines/split-pipe/split-pipe_mkref/main'
include { PARSEBIO_PIPELINE } from '../modules/pipelines/split-pipe/split-pipe_all/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR PARSE BIOSCIENCES DATA
        Based on the sample wells, the fastq sequences must be split 
        based on the first barcode round. Performs both custom- and 
        commercial pre-processing (using Parse Biosciences pipeline v1.3.1) 
        enabling comparison between the methods and validation of steps. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow parse_workflow {
    take:
        ch_samplesheet
    main:       
        // Demultiplex the fastq files based on the sample wells
        PARSEBIO_PIPELINE_DEMUX(ch_samplesheet)
        
        // Only run Parse pipeline if the path is defined and exists
        if (params.splitpipe_installation && file(params.splitpipe_installation).exists()) {
            PARSEBIO_PIPELINE_MKREF()
            PARSEBIO_PIPELINE(PARSEBIO_PIPELINE_DEMUX.out.splitted_files, PARSEBIO_PIPELINE_MKREF.out)
        } else {
            log.warn "Parse Biosciences pipeline directory not provided or doesn't exist: '${params.splitpipe_installation}'"
        }
    
    emit:
        // Result: demultiplexed fastq files from split-pipe
        PARSEBIO_PIPELINE_DEMUX.out.splitted_files
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
