//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RM_VARBASES } from '../../../modules/local/tools/cutadapt/main'
include { BDRHAP_PIPELINE_MKREF } from '../../../modules/local/pipelines/rhapsody_pipeline/rhapsody_mkref/main'
include { BDRHAP_PIPELINE_YAML } from '../../../modules/local/pipelines/rhapsody_pipeline/rhapsody_create_yml/main'
include { BDRHAP_PIPELINE } from '../../../modules/local/pipelines/rhapsody_pipeline/rhapsody_full/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR BD RHAPSODY DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow bd_rhapsody_workflow {
    take:
        ch_samplesheet
    main:
        // Remove variable bases (0-3) from the fastq files using cutadapt
        RM_VARBASES(ch_samplesheet)

        // Only run BD Rhapsody pipeline if the path is defined and exists
        if (params.rhapsody_installation) {
            BDRHAP_PIPELINE_MKREF()

            // Use .first() to treat the reference as a reusable Value Channel
            BDRHAP_PIPELINE_YAML(ch_samplesheet, BDRHAP_PIPELINE_MKREF.out.reference.first())

            // Join all 4 outputs from YAML generation
            BDRHAP_PIPELINE_YAML.out.run_name
                .join(BDRHAP_PIPELINE_YAML.out.bd_ref_path)
                .join(BDRHAP_PIPELINE_YAML.out.yaml_file)
                .join(BDRHAP_PIPELINE_YAML.out.samplesheet)
                .multiMap { meta, run_name, ref, yaml, sheet ->
                    run_name_ch: [meta, run_name]
                    ref_ch:      [meta, ref]
                    yaml_ch:     [meta, yaml]
                    sheet_ch:    [meta, sheet]
                }
                .set { ch_bd_pipeline_inputs }

            BDRHAP_PIPELINE(
                ch_bd_pipeline_inputs.run_name_ch,
                ch_bd_pipeline_inputs.ref_ch,
                ch_bd_pipeline_inputs.yaml_ch,
                ch_bd_pipeline_inputs.sheet_ch
            )

        } else {
            log.warn "BD Rhapsody pipeline directory not provided or doesn't exist: '${params.rhapsody_installation}'"
        }

    emit:
        data_output = RM_VARBASES.out.trimmed_files
        bc_whitelist    = params.bc_whitelist
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
