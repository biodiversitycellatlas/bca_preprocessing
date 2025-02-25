// ===============  Parse Biosciences  ==============  \\ 
// Based on the sample wells, the fastq sequences      \\
// must be split based on the first barcode round.     \\
// Performs both custom- and commercial pre-processing \\
// (using the Parse Biosciences pipeline v1.3.1)       \\
// enabling comparison between the methods and         \\
// validation of steps.       

include { DOWNLOAD_DATA } from '../modules/download'

include { DEMUX_SPIPE } from '../modules/demultiplex'
include { DEMUX_UMITOOLS } from '../modules/demux_umitools'

include { REFGEN_PARSEBIO } from '../modules/refgen_parsebio'
include { MAPPING_PARSEBIO } from '../modules/mapping_parsebio'

workflow parse_workflow {
    take:
        sample_ids
    main:
        DOWNLOAD_DATA(sample_ids)
        
        groups = Channel.fromList(params.groups)
        comb_data = DOWNLOAD_DATA.out.fastq_files.combine(groups)

        DEMUX_SPIPE(comb_data)
        DEMUX_UMITOOLS(comb_data)
        
        REFGEN_PARSEBIO()
        MAPPING_PARSEBIO(DEMUX_SPIPE.out.splitted_files, REFGEN_PARSEBIO.out)
    emit:
        DEMUX_UMITOOLS.out.splitted_files
}
