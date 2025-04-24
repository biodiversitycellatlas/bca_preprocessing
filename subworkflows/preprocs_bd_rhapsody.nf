//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_DATA } from '../modules/download'
include { RM_VARBASES } from '../modules/rm_varbases_bdrhap'
include { BDRHAP_PIPELINE } from '../modules/bdrhap_pipeline'
include { BDRHAP_PIPELINE_MKREF } from '../modules/bdrhap_pipeline_mkref'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR BD RHAPSODY DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow bd_rhapsody_workflow {
    take:
        ch_samplesheet
    main:
        // Import the fastq files into the nf workdir using sym links to the original files
        DOWNLOAD_DATA(ch_samplesheet)     

        // Remove variable bases (0-3) from the fastq files using cutadapt
        RM_VARBASES(DOWNLOAD_DATA.out)

        // Only run BD Rhapsody pipeline if the path is defined and exists
        if (params.bdrhap_pipeline_dir && file(params.bdrhap_pipeline_dir).exists()) {
            BDRHAP_PIPELINE_MKREF()
            BDRHAP_PIPELINE(ch_samplesheet, BDRHAP_PIPELINE_MKREF.out)
        } else {
            log.warn "BD Rhapsody pipeline directory not provided or doesn't exist: '${params.bdrhap_pipeline_dir}'"
        }

    emit:
        RM_VARBASES.out
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/