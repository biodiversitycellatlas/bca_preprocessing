include { DOWNLOAD_DATA } from '../modules/download'
include { CR_PIPELINE_MKREF } from '../modules/cellranger_pipeline_mkref'
include { CR_PIPELINE } from '../modules/cellranger_pipeline'

workflow oak_seq_workflow {
    take:
        sample_ids
    main:
        // Import the fastq files into the nf workdir using sym links to the original files
        DOWNLOAD_DATA(sample_ids)

        // Create reference for CellRanger pipeline
        CR_PIPELINE_MKREF()

        // Run CellRanger pipeline
        CR_PIPELINE(CR_PIPELINE_MKREF.out)

    emit:
        DOWNLOAD_DATA.out
}
