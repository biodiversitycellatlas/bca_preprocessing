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
            BDRHAP_PIPELINE_YAML(ch_samplesheet, BDRHAP_PIPELINE_MKREF.out)
            BDRHAP_PIPELINE(BDRHAP_PIPELINE_YAML.out.run_name, BDRHAP_PIPELINE_YAML.out.bd_ref_path, BDRHAP_PIPELINE_YAML.out.yaml_file)
        } else {
            log.warn "BD Rhapsody pipeline directory not provided or doesn't exist: '${params.rhapsody_installation}'"
        }

    emit:
        RM_VARBASES.out
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
