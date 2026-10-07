//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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
        def ch_versions = Channel.empty()
        def ch_cellranger_outs = Channel.empty()

        // Only run Cell Ranger pipeline if perform_cellranger is set to true
        if (params.perform_cellranger) {
            CR_PIPELINE_MKREF()

            // Use .first() to allow the reference index to be reused for all samples
            CR_PIPELINE(ch_samplesheet, CR_PIPELINE_MKREF.out.reference.first())
            ch_cellranger_outs = CR_PIPELINE.out.outs
            ch_versions = ch_versions.mix(CR_PIPELINE_MKREF.out.versions, CR_PIPELINE.out.versions)
        }

    emit:
        data_output     = ch_samplesheet
        cellranger_outs = ch_cellranger_outs
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
