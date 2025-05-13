//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { mapping_starsolo_workflow } from '../subworkflows/mapping_starsolo'
include { mapping_alevin_workflow   } from '../subworkflows/mapping_alevin'

include { FASTQC                    } from '../modules/fastqc'
include { KRAKEN_CREATE_DB          } from '../modules/kraken_create_db'
include { KRAKEN                    } from '../modules/kraken'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow QC_mapping_workflow {
    take:
        data_output
        bc_whitelist

    main:
        // Quality Control
        FASTQC(data_output)

        // Combine outputs from all processes into a single channel
        all_outputs = FASTQC.out

        // Mapping: STARsolo, Alevin-fry, or both
        if (params.mapping_software == "starsolo") {
            mapping_starsolo_workflow(data_output, bc_whitelist, all_outputs)
            mapping_files = mapping_starsolo_workflow.out.mapping_files
            all_outputs = all_outputs.mix(mapping_starsolo_workflow.out.all_outputs)

        } else if (params.mapping_software == "alevin") {
            mapping_alevin_workflow(data_output, bc_whitelist, all_outputs)
            mapping_files = mapping_alevin_workflow.out.mapping_files
            all_outputs = all_outputs.mix(mapping_alevin_workflow.out.all_outputs)

        } else if (params.mapping_software == "both") {
            mapping_starsolo_workflow(data_output, bc_whitelist, all_outputs)
            mapping_alevin_workflow(data_output, bc_whitelist, all_outputs)
            mapping_files = mapping_alevin_workflow.out.mapping_files.mix(mapping_starsolo_workflow.out.mapping_files)
            all_outputs = all_outputs.mix(mapping_alevin_workflow.out.all_outputs)
            all_outputs = all_outputs.mix(mapping_starsolo_workflow.out.all_outputs)
        
        } else {
            error "Invalid mapping software specified. Use one of the following parameters: 'starsolo', 'alevin' or 'both'."
        }

        // Inspecting unmapped reads
        // Conditionally run Kraken only if params.perform_kraken is true and DB path is valid
        if (params.perform_kraken && file(params.kraken_db_path).exists()) {
            KRAKEN_CREATE_DB()
            KRAKEN(KRAKEN_CREATE_DB.out.db_path_file, mapping_files) 
            all_outputs = all_outputs.mix(KRAKEN.out)
        } else {
            log.info "Skipping Kraken steps as 'perform_kraken' is false."
        }

    emit:
        mapping_files
        all_outputs 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/