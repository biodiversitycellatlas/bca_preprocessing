//
// Subworkflow with functionality specific to the workflow 'mapping_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SALMON_SPLICI     } from '../../../modules/local/tools/salmon_alevin/salmon_splici/main'
include { SALMON_INDEX      } from '../../../modules/local/tools/salmon_alevin/salmon_index/main'
include { ALEVIN_FRY        } from '../../../modules/local/tools/salmon_alevin/alevin-fry/main'
include { ALEVIN_QC         } from '../../../modules/local/tools/alevinQC/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN ALEVIN-FRY MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow mapping_alevin_workflow {
    take:
        data_output
        bc_whitelist

    main:
        // Build one Salmon splici and index for all samples
        SALMON_SPLICI(data_output.first())
        SALMON_INDEX(SALMON_SPLICI.out.splici_index_reference)

        // Run Alevin-fry for each sample
        ALEVIN_FRY(
            data_output,
            bc_whitelist,
            SALMON_SPLICI.out.splici_index_reference,
            SALMON_INDEX.out.salmon_index
        )

        // Quality control report for Alevin-fry outputs
        ALEVIN_QC(ALEVIN_FRY.out.mapping_files)

    emit:
        mapping_files   = ALEVIN_FRY.out.mapping_files
        af_meta_info    = ALEVIN_FRY.out.af_meta_info
        af_quant_json   = ALEVIN_FRY.out.af_quant_json
        af_cell_meta    = ALEVIN_FRY.out.af_cell_meta
        qc_reports      = ALEVIN_QC.out.alevinQC_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
