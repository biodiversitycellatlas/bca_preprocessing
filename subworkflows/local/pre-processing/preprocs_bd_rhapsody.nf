//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RM_VARBASES } from '../../../modules/local/tools/cutadapt/main'
include { BDRHAP_PIPELINE } from '../../../modules/local/pipelines/rhapsody_pipeline/rhapsody_full/main'
include { BDRHAP_PIPELINE_MKREF } from '../../../modules/local/pipelines/rhapsody_pipeline/rhapsody_mkref/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR BD RHAPSODY DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow bd_rhapsody_workflow {
    take:
        ch_samplesheet
    main:
        // Remove variable bases (0-3) from the fastq files using cutadapt
        RM_VARBASES(ch_samplesheet)

        // Only run BD Rhapsody pipeline if the path is defined and exists
        if (params.rhapsody_installation && file(params.rhapsody_installation).exists()) {
            BDRHAP_PIPELINE_MKREF()
            BDRHAP_PIPELINE(ch_samplesheet, BDRHAP_PIPELINE_MKREF.out)
        } else {
            log.warn "BD Rhapsody pipeline directory not provided or doesn't exist: '${params.rhapsody_installation}'"
        }

    emit:
        RM_VARBASES.out
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/