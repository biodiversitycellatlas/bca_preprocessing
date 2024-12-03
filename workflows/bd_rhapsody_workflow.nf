// ===================  BD Rhapsody ================  \\ 
// A basic approach, where first quality control is   \\
// performed and mapping using the (complete) fastq   \\
// sequences.                                         \\

include { DOWNLOAD_DATA } from '../modules/download'

workflow bd_rhapsody_workflow {
    take:
        sample_ids
    main:
        DOWNLOAD_DATA(sample_ids)        
    emit:
        DOWNLOAD_DATA.out
}
