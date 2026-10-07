//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MARSSEQ_BUILD_READS } from '../../../modules/local/custom/demultiplex/marsseq_demultiplex/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR MARS-SEQ DATA (v1 & v2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow marsseq_workflow {
    take:
        ch_samplesheet

    main:
        // Define channels
        out_samplesheet = Channel.empty()
        ch_versions = Channel.empty()

        if (params.perform_demultiplexing) {
            log.info "MARS-seq '${params.protocol}': rebuilding reads, the batch barcode moves in front of the cell barcode and UMI."
            MARSSEQ_BUILD_READS(ch_samplesheet)
            out_samplesheet = MARSSEQ_BUILD_READS.out.reformatted_files
            ch_versions = ch_versions.mix(MARSSEQ_BUILD_READS.out.versions)
        } else {
            log.info "Skipping demultiplexing as perform_demultiplexing is set to false"
            out_samplesheet = ch_samplesheet
        }

    emit:
        data_output = out_samplesheet
        versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
