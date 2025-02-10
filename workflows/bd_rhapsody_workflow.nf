// ===================  BD Rhapsody ================  \\ 
// A basic approach, where first quality control is   \\
// performed and mapping using the (complete) fastq   \\
// sequences.                                         \\

include { DOWNLOAD_DATA } from '../modules/download'
include { RM_VARBASES } from '../modules/rm_varbases_bdrhap'

workflow bd_rhapsody_workflow {
    take:
        sample_ids
    main:
        DOWNLOAD_DATA(sample_ids)     
        RM_VARBASES(DOWNLOAD_DATA.out)   
    emit:
        RM_VARBASES.out.fastq_noVB_files
}
