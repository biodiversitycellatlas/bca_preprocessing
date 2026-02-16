//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { mapping_starsolo_workflow } from '../subworkflows/local/mapping/mapping_starsolo'
include { mapping_alevin_workflow   } from '../subworkflows/local/mapping/mapping_alevin'

include { FASTQC                    } from '../modules/local/tools/fastqc/main'
include { KRAKEN_CREATE_DB          } from '../modules/local/tools/kraken/kraken_create_db/main'
include { KRAKEN                    } from '../modules/local/tools/kraken/kraken_classify/main'
include { KRONA                     } from '../modules/local/tools/krona/main'
include { PAVIAN                    } from '../modules/local/tools/pavian/main'


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

        // Initialize variables to ensure scope visibility
        def ch_all_outputs   = FASTQC.out
        def ch_mapping_files = Channel.empty()
        def ch_starsolo_bam  = Channel.empty()
        def ch_starsolo_genefull50_raw  = Channel.empty()

        // Mapping: STARsolo, Alevin-fry, or both
        if (params.mapping_software == "starsolo") {
            mapping_starsolo_workflow(data_output, bc_whitelist, ch_all_outputs)

            ch_mapping_files         =  mapping_starsolo_workflow.out.mapping_files
            ch_all_outputs           =  ch_all_outputs.mix(mapping_starsolo_workflow.out.all_outputs)
            ch_starsolo_bam          =  mapping_starsolo_workflow.out.starsolo_bam
            ch_starsolo_genefull50_raw  = mapping_starsolo_workflow.out.starsolo_genefull50_raw

        } else if (params.mapping_software == "alevin") {
            mapping_alevin_workflow(data_output, bc_whitelist, ch_all_outputs)
            ch_mapping_files = mapping_alevin_workflow.out.mapping_files
            ch_all_outputs = ch_all_outputs.mix(mapping_alevin_workflow.out.all_outputs)

        } else if (params.mapping_software == "both") {
            mapping_starsolo_workflow(data_output, bc_whitelist, ch_all_outputs)
            mapping_alevin_workflow(data_output, bc_whitelist, ch_all_outputs)
            ch_mapping_files = mapping_alevin_workflow.out.mapping_files.mix(mapping_starsolo_workflow.out.mapping_files)
            ch_all_outputs = ch_all_outputs.mix(mapping_alevin_workflow.out.all_outputs)
            ch_all_outputs = ch_all_outputs.mix(mapping_starsolo_workflow.out.all_outputs)
            ch_starsolo_bam  = mapping_starsolo_workflow.out.starsolo_bam
            ch_starsolo_genefull50_raw  = mapping_starsolo_workflow.out.starsolo_genefull50_raw

        } else {
            error "Invalid mapping software specified. Use one of the following parameters: 'starsolo', 'alevin' or 'both'."
        }

        // Inspecting unmapped reads
        // Conditionally run Kraken only if params.perform_kraken is true
        if (params.perform_kraken && params.mapping_software != "alevin") {
            KRAKEN_CREATE_DB()
            KRAKEN(KRAKEN_CREATE_DB.out.db_path_file, ch_starsolo_bam)
            ch_all_outputs = ch_all_outputs.mix(KRAKEN.out)
            KRONA(KRAKEN.out)
            PAVIAN(KRAKEN.out)
        }

    emit:
        mapping_files = ch_mapping_files
        all_outputs   = ch_all_outputs
        starsolo_genefull50_raw  = ch_starsolo_genefull50_raw
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
