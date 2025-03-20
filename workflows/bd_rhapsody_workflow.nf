include { DOWNLOAD_DATA } from '../modules/download'
include { RM_VARBASES } from '../modules/rm_varbases_bdrhap'
include { DEMUX_UMITOOLS_BDRHAP } from '../modules/demux_umitools_bdrhap'
include { CONV_3CB_INDEX } from '../modules/conv_3cb_index'

include { BDRHAP_PIPELINE } from '../modules/bdrhap_pipeline'
include { BDRHAP_PIPELINE_MKREF } from '../modules/bdrhap_pipeline_mkref'


workflow bd_rhapsody_workflow {
    take:
        sample_ids
    main:
        // Import the fastq files into the nf workdir using sym links to the original files
        DOWNLOAD_DATA(sample_ids)     

        // Remove variable bases (0-3) from the fastq files using cutadapt
        RM_VARBASES(DOWNLOAD_DATA.out)

        // Extract the CBs and UMI's and add them to the header of the cDNA file
        // DEMUX_UMITOOLS_BDRHAP(RM_VARBASES.out)

        // Convert 3CBs to an index, to compare cells with the BD rhapsody pipeline
        // CONV_3CB_INDEX(DEMUX_UMITOOLS_BDRHAP.out.splitted_files)

        // Create the reference files for the BD Rhapsody pipeline
        BDRHAP_PIPELINE_MKREF()

        // Run the BD Rhapsody pipeline
        BDRHAP_PIPELINE(BDRHAP_PIPELINE_MKREF.out)

    emit:
        RM_VARBASES.out
}
