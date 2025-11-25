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
include { ALEVIN_QC          } from '../../../modules/local/tools/alevinQC/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN ALEVIN-FRY MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow mapping_alevin_workflow {
    take:
        data_output
        bc_whitelist
        all_outputs

    main:
        SALMON_SPLICI(data_output)
        SALMON_INDEX(data_output, SALMON_SPLICI.out)
        ALEVIN_FRY(data_output, bc_whitelist, SALMON_SPLICI.out, SALMON_INDEX.out)

        ALEVIN_QC(data_output, ALEVIN_FRY.out)

    emit:
        mapping_files = ALEVIN_FRY.out
        all_outputs
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
