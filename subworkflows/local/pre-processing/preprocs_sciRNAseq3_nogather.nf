/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SCIROCKET_DEMUX       } from '../../../modules/local/tools/scirocket/scirocket_demux/main'
include { FASTP                 } from '../../../modules/local/tools/fastp/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR SCI-RNA-SEQ3 DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow sciRNAseq3_nogather_workflow {
    take:
        ch_samplesheet

    main:
        // Define channels
        out_whitelist = Channel.empty()
        out_samplesheet = Channel.empty()

        if (params.perform_demultiplexing) {
            log.info "Starting demultiplexing with sci-rocket"
            SCIROCKET_DEMUX(ch_samplesheet)
            out_samplesheet = SCIROCKET_DEMUX.out.demux_samplesheet
            out_whitelist = SCIROCKET_DEMUX.out.bc_whitelist
        } else {
            log.info "Skipping demultiplexing as perform_demultiplexing is set to false"
            out_samplesheet = ch_samplesheet
        }

        // Trimming adapters and low-quality reads
        FASTP(out_samplesheet)

    emit:
        data_output     = FASTP.out.trimmed_files
        bc_whitelist    = out_whitelist
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
