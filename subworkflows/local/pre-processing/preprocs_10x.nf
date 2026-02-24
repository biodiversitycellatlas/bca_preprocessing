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

        // Convert row-map -> tuple(meta, R2, R1, indices, samplesheet_file)
        ch_samplesheet_tuple = ch_samplesheet.map { row ->
            tuple(
                row.meta,
                file(row.fastq_cDNA),
                file(row.fastq_BC_UMI),
                row.fastq_indices ? file(row.fastq_indices) : [],
                file(row.input_file)
            )
        }

        // Only run Cell Ranger pipeline if perform_cellranger is set to true
        if (params.perform_cellranger) {
            CR_PIPELINE_MKREF()

            // Use .first() to allow the reference index to be reused for all samples
            CR_PIPELINE(ch_samplesheet_tuple, CR_PIPELINE_MKREF.out.first())
        }

    emit:
        data_output     = ch_samplesheet_tuple
        bc_whitelist    = DOWNLOAD_DATA.out
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
