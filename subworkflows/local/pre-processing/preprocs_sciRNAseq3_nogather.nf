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

        // The mappers' whitelists are written by sci-rocket, one set per sample, and the
        // barcode read they describe is rebuilt by it: without demultiplexing there is neither
        if (!params.perform_demultiplexing) {
            error("protocol 'sciRNAseq3' needs 'perform_demultiplexing = true': sci-rocket rebuilds the 48 nt barcode read that the mappers read, and writes each sample's barcode whitelists.")
        }
        log.info "Starting demultiplexing with sci-rocket"
        if (!params.bc_whitelist) {
            error("sci-rocket demultiplexes on the barcodes in 'bc_whitelist', which is not set.")
        }
        SCIROCKET_DEMUX(ch_samplesheet, bc_whitelist)
        out_samplesheet = SCIROCKET_DEMUX.out.demux_samplesheet

        // Each sample is mapped against its own p5, p7 and RT barcodes (sci-rocket itself
        // gives every sample the union of the run's): the reads hold only those after
        // demultiplexing, so the result is the same, and STARsolo, which expands the segment
        // lists into every combination, stays smaller
        out_whitelist = SCIROCKET_DEMUX.out.bc_whitelist
            .map { meta, p7, p5, ligation, rt -> [meta.id, [p7, p5, ligation, rt]] }
        ch_versions = ch_versions.mix(SCIROCKET_DEMUX.out.versions)

        // perform_trimming if set, otherwise the protocol's default (conf/seqtech_parameters.config)
        def perform_trimming = params.perform_trimming != null
            ? params.perform_trimming
            : params.seqtech_parameters[params.protocol].perform_trimming
        if (perform_trimming) {
            // Trimming adapters and low-quality reads
            FASTP(out_samplesheet)
            out_data = FASTP.out.trimmed_files
            ch_versions = ch_versions.mix(FASTP.out.versions)
        } else {
            log.info "Skipping trimming: perform_trimming is false${params.perform_trimming == null ? " (the protocol's default)" : ""}"
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
