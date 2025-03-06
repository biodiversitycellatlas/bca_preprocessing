// ===================  OAK seq ================  \\ 
// A basic approach, where first quality control is   \\
// performed and mapping using the (complete) fastq   \\
// sequences.                                         \\

include { DOWNLOAD_DATA } from '../modules/download'

workflow oak_seq_workflow {
    take:
        sample_ids
    main:
        // Import the fastq files into the nf workdir using sym links to the original files
        DOWNLOAD_DATA(sample_ids)     

    emit:
        DOWNLOAD_DATA.out
}
