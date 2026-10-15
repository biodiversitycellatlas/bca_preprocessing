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
        bc_whitelist

    main:
        // Define channels
        out_whitelist = channel.empty()
        out_samplesheet = channel.empty()
        out_data = channel.empty()
        ch_versions = channel.empty()

        if (params.perform_demultiplexing) {
            log.info "Starting demultiplexing with sci-rocket"
            if (!params.bc_whitelist) {
                error("sci-rocket demultiplexes on the barcodes in 'bc_whitelist', which is not set.")
            }
            SCIROCKET_DEMUX(ch_samplesheet, bc_whitelist)
            out_samplesheet = SCIROCKET_DEMUX.out.demux_samplesheet
            out_whitelist = SCIROCKET_DEMUX.out.bc_whitelist
            ch_versions = ch_versions.mix(SCIROCKET_DEMUX.out.versions)
        } else {
            log.info "Skipping demultiplexing as perform_demultiplexing is set to false"
            out_samplesheet = ch_samplesheet
        }

        if (params.perform_trimming) {
            // Trimming adapters and low-quality reads
            FASTP(out_samplesheet)
            out_data = FASTP.out.trimmed_files
            ch_versions = ch_versions.mix(FASTP.out.versions)
        } else {
            log.info "Skipping trimming as perform_trimming is set to false"
            out_data = out_samplesheet
        }

    emit:
        data_output     = out_data
        bc_whitelist    = out_whitelist
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
