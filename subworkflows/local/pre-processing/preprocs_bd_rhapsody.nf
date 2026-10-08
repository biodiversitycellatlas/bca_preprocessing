//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RM_VARBASES } from '../../../modules/local/tools/cutadapt/main'
include { STARSOLO_INDEX as BDRHAP_STAR_INDEX } from '../../../modules/local/tools/star/starsolo_genome_generate/main'
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
        def ch_versions = channel.empty()

        // Enhanced Beads reads start with 0-3 variable bases, the V1 beads have a fixed layout
        if (params.protocol == 'bd_rhapsody_enhancedbeads') {
            // Remove variable bases (0-3) from the fastq files using cutadapt
            RM_VARBASES(ch_samplesheet)
            preprocs_samplesheet = RM_VARBASES.out.trimmed_files
            ch_versions = ch_versions.mix(RM_VARBASES.out.versions)

        } else {
            log.info "Skipping variable base removal: '${params.protocol}' reads do not carry variable bases."
            preprocs_samplesheet = ch_samplesheet
        }

        // Only run BD Rhapsody pipeline if the path is defined and exists
        if (params.rhapsody_installation && file(params.rhapsody_installation).exists()) {
            // One reference for all samples, its sjdbOverhang set from the longest cDNA read
            // Assigned without 'def': Nextflow 24.04 fails to compile it here otherwise
            all_cdna_ch = ch_samplesheet
                .map { meta, fastq_cDNA, fastq_BC_UMI, fastq_indices, input_file -> fastq_cDNA }
                .collect()
            BDRHAP_STAR_INDEX(all_cdna_ch, file(params.ref_gtf), file(params.ref_fasta), '_bdrhapsody')

            // Packs the index and GTF into the archive the BD Rhapsody pipeline reads
            BDRHAP_PIPELINE_MKREF(BDRHAP_STAR_INDEX.out.index, file(params.ref_gtf))

            // The BD Rhapsody pipeline handles the bead layout itself, so it gets the untrimmed reads
            // Use .first() to treat the reference as a reusable Value Channel
            BDRHAP_PIPELINE_YAML(ch_samplesheet, BDRHAP_PIPELINE_MKREF.out.reference.first())

            // Reference, yaml and fastq files are emitted together, keeping them in sync per sample
            BDRHAP_PIPELINE(BDRHAP_PIPELINE_YAML.out.pipeline_input)

            ch_versions = ch_versions.mix(
                BDRHAP_STAR_INDEX.out.versions,
                BDRHAP_PIPELINE_MKREF.out.versions,
                BDRHAP_PIPELINE_YAML.out.versions,
                BDRHAP_PIPELINE.out.versions
            )

        } else {
            log.warn "BD Rhapsody pipeline directory not provided or doesn't exist: '${params.rhapsody_installation}'"
        }

    emit:
        data_output = preprocs_samplesheet
        versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
