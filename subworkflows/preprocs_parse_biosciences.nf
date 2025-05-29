//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_DATA } from '../modules/download'
include { PARSEBIO_CUSTOM_DEMUX } from '../modules/parsebio_custom_demux'
include { PARSEBIO_PIPELINE_DEMUX } from '../modules/parsebio_pipeline_demux'
include { PARSEBIO_PIPELINE_MKREF } from '../modules/parsebio_pipeline_mkref'
include { PARSEBIO_PIPELINE } from '../modules/parsebio_pipeline'

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
        PARSEBIO_CUSTOM_DEMUX(ch_samplesheet)
        
        // Only run Parse pipeline if the path is defined and exists
        if (params.parsebio_pipeline_dir && file(params.parsebio_pipeline_dir).exists()) {
            PARSEBIO_PIPELINE_MKREF()
            PARSEBIO_PIPELINE(PARSEBIO_PIPELINE_DEMUX.out.splitted_files, PARSEBIO_PIPELINE_MKREF.out)
        } else {
            log.warn "Parse Biosciences pipeline directory not provided or doesn't exist: '${params.parsebio_pipeline_dir}'"
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
