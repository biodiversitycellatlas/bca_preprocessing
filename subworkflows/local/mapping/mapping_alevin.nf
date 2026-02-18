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
        SALMON_SPLICI(data_output)
        SALMON_INDEX(SALMON_SPLICI.out.splici_index_reference)

        // Join SALMON_INDEX outputs with SALMON_SPLICI outputs
        SALMON_INDEX.out.salmon_index
            .join(data_output)
            .join(SALMON_SPLICI.out.splici_index_reference)
            .multiMap { meta, index, cdna, bc_umi, indices, input_file, splici_ref ->
                ssheet_ch:      [meta, cdna, bc_umi, indices, input_file]
                index_ch:       [meta, index]
                splici_ref_ch:  [meta, splici_ref]
            }
            .set { ch_alevin_inputs }

        ALEVIN_FRY(
            ch_alevin_inputs.ssheet_ch,
            bc_whitelist,
            ch_alevin_inputs.splici_ref_ch,
            ch_alevin_inputs.index_ch
        )

        // Quality control report for Alevin-fry outputs
        ALEVIN_QC(ALEVIN_FRY.out.mapping_files)

    emit:
        mapping_files   = ALEVIN_FRY.out.mapping_files
        qc_reports      = ALEVIN_QC.out.alevinQC_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
