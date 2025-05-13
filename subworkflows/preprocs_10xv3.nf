//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_DATA } from '../modules/download'
include { CR_PIPELINE_MKREF } from '../modules/cellranger_pipeline_mkref'
include { CR_PIPELINE } from '../modules/cellranger_pipeline'

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
        if (params.cellranger_dir && file(params.cellranger_dir).exists()) {
            CR_PIPELINE_MKREF()
            CR_PIPELINE(ch_samplesheet, CR_PIPELINE_MKREF.out)
        } else {
            log.warn "Cell Ranger pipeline directory not provided or doesn't exist: '${params.cellranger_dir}'"
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