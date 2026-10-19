//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTP                 } from '../../../modules/local/tools/fastp/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING WHEN GIVEN A SEQSPEC FILE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow seqspec_workflow {
    take:
        ch_samplesheet
    main:
        // Define channels
        out_data = channel.empty()
        ch_versions = channel.empty()

        if (params.perform_trimming) {
            // Trimming adapters and low-quality reads
            FASTP(ch_samplesheet)
            out_data = FASTP.out.trimmed_files
            ch_versions = ch_versions.mix(FASTP.out.versions)
        } else {
            log.info "Skipping trimming as perform_trimming is set to false"
            out_data = ch_samplesheet
        }

    emit:
        data_output     = out_data
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
