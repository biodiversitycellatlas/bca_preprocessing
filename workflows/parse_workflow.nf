// ===============  Parse Biosciences  ==============  \\ 
// Based on the sample wells, the fastq sequences      \\
// must be split based on the first barcode round.     \\
// Performs both custom- and commercial pre-processing \\
// (using the Parse Biosciences pipeline v1.3.1)       \\
// enabling comparison between the methods and         \\
// validation of steps.       

include { DOWNLOAD_DATA } from '../modules/download'
include { DEMULTIPLEX } from '../modules/demultiplex'
include { REFGEN_PARSEBIO } from '../modules/refgen_parsebio'
include { MAPPING_PARSEBIO } from '../modules/mapping_parsebio'

workflow parse_workflow {
    take:
        sample_ids
    main:
        DOWNLOAD_DATA(sample_ids)
        DEMULTIPLEX(DOWNLOAD_DATA.out)
        
        REFGEN_PARSEBIO()
        MAPPING_PARSEBIO(DEMULTIPLEX.out.splitted_files, REFGEN_PARSEBIO.out)
    emit:
        DEMULTIPLEX.out.splitted_files
}
