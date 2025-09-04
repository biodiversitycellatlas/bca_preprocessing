//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_DATA } from '../../../modules/local/custom/manipulate/download_files/main'
include { CR_PIPELINE_MKREF } from '../../../modules/local/pipelines/cellranger/cellranger_mkref/main'
include { CR_PIPELINE } from '../../../modules/local/pipelines/cellranger/cellranger_count/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR 10X GENOMICS DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow tenx_genomics_workflow {
    take:
        ch_samplesheet
    main:
        // Download data from the specified path, in this case the barcode whitelist
        DOWNLOAD_DATA()

        // Only run Cell Ranger pipeline if the path is defined and exists
        if (params.perform_cellranger) {
            CR_PIPELINE_MKREF()
            CR_PIPELINE(ch_samplesheet, CR_PIPELINE_MKREF.out)
        }

    emit:
        data_output     = ch_samplesheet
        bc_whitelist    = DOWNLOAD_DATA.out
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
