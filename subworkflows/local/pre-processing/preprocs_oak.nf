//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { OAK_DEMUX         } from '../../../modules/local/custom/demultiplex/oak_demultiplex/main'
include { CR_PIPELINE_MKREF } from '../../../modules/local/pipelines/cellranger/cellranger_mkref/main'
include { CR_PIPELINE       } from '../../../modules/local/pipelines/cellranger/cellranger_count/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR OAK (v1) DATA
        Every samplesheet row is one OAK aliquot, demultiplexed from the shared
        Undetermined FASTQs on its i7 (optional) and i5 index. Each aliquot is then
        processed as a 10x 3' v3.1 sample.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow oak_workflow {
    take:
        ch_samplesheet

    main:
        // Define channels
        out_samplesheet = Channel.empty()
        ch_versions = Channel.empty()
        ch_cellranger_outs = Channel.empty()

        if (params.perform_demultiplexing) {
            log.info "OAK '${params.protocol}': demultiplexing aliquots on their i7 (optional) and i5 index."
            OAK_DEMUX(ch_samplesheet)
            out_samplesheet = OAK_DEMUX.out.demux_files
            ch_versions = ch_versions.mix(OAK_DEMUX.out.versions)
        } else {
            log.info "Skipping demultiplexing as perform_demultiplexing is set to false"
            out_samplesheet = ch_samplesheet
        }

        // Only run Cell Ranger pipeline if perform_cellranger is set to true
        if (params.perform_cellranger) {
            CR_PIPELINE_MKREF()

            // Use .first() to allow the reference index to be reused for all samples
            CR_PIPELINE(out_samplesheet, CR_PIPELINE_MKREF.out.reference.first())
            ch_cellranger_outs = CR_PIPELINE.out.outs
            ch_versions = ch_versions.mix(CR_PIPELINE_MKREF.out.versions, CR_PIPELINE.out.versions)
        }

    emit:
        data_output     = out_samplesheet
        cellranger_outs = ch_cellranger_outs
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
