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
        // Enhanced Beads reads start with 0-3 variable bases, the V1 beads have a fixed layout
        if (params.protocol == 'bd_rhapsody_enhancedbeads') {
            // Remove variable bases (0-3) from the fastq files using cutadapt
            RM_VARBASES(ch_samplesheet)
            preprocs_samplesheet = RM_VARBASES.out.trimmed_files

        } else {
            log.info "Skipping variable base removal: '${params.protocol}' reads do not carry variable bases."
            preprocs_samplesheet = ch_samplesheet
        }

        // Only run BD Rhapsody pipeline if the path is defined and exists
        if (params.rhapsody_installation && file(params.rhapsody_installation).exists()) {
            BDRHAP_PIPELINE_MKREF()

            // The BD Rhapsody pipeline handles the bead layout itself, so it gets the untrimmed reads
            // Use .first() to treat the reference as a reusable Value Channel
            BDRHAP_PIPELINE_YAML(ch_samplesheet, BDRHAP_PIPELINE_MKREF.out.reference.first())

            // Reference, yaml and fastq files are emitted together, keeping them in sync per sample
            BDRHAP_PIPELINE(BDRHAP_PIPELINE_YAML.out.pipeline_input)

        } else {
            log.warn "BD Rhapsody pipeline directory not provided or doesn't exist: '${params.rhapsody_installation}'"
        }

    emit:
        data_output = preprocs_samplesheet
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
