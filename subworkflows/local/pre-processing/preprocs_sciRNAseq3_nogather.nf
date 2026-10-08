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
        out_whitelist = channel.empty()
        out_samplesheet = channel.empty()
        ch_versions = channel.empty()

        if (params.perform_demultiplexing) {
            log.info "Starting demultiplexing with sci-rocket"
            SCIROCKET_DEMUX(ch_samplesheet)
            out_samplesheet = SCIROCKET_DEMUX.out.demux_samplesheet
            out_whitelist = SCIROCKET_DEMUX.out.bc_whitelist
            ch_versions = ch_versions.mix(SCIROCKET_DEMUX.out.versions)
        } else {
            log.info "Skipping demultiplexing as perform_demultiplexing is set to false"
            out_samplesheet = ch_samplesheet
        }

        // Trimming adapters and low-quality reads
        FASTP(out_samplesheet)
        ch_versions = ch_versions.mix(FASTP.out.versions)

    emit:
        data_output     = FASTP.out.trimmed_files
        bc_whitelist    = out_whitelist
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
