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
        data = DOWNLOAD_DATA(sample_ids)
        
        // include BCA00.. in both channels to not have to split nvec. Match on these then and it should work
        // Nvec - BCA001/2
        groups = Channel.of(['ACMEsorb_GM', 'A1-A3'], ['EMA_ACME_GM', 'A4-A6'], ['DSP', 'A7-A9'], ['Parse_fix', 'A10-A12'])
        // Nvec - BCA003/4
        // groups = Channel.of(['ACMEsorb_cold_PR', 'A1-A3'], ['Vivophix_sonic', 'A4-A6'], ['DSP_CMFSW', 'A7-A9'])
        // Nvec nuclei - BCA007/8
        // groups = Channel.of(['DCE_PRS', 'A1-A3'], ['DCE_FH', 'A10-A12'])

        // Tcas - BCA003/4
        // groups = Channel.of(['all', 'A10-A12'], ['A10', 'A10'], ['A11', 'A11'], ['A12', 'A12'])
        // Tcas nuclei - BCA007/8
        // groups = Channel.of(['DCE_PRS', 'A4-A6'], ['DCE_FH', 'A7-A9'])

        comb_data = data.combine(groups)
        DEMULTIPLEX(comb_data)
        
        REFGEN_PARSEBIO()
        MAPPING_PARSEBIO(DEMULTIPLEX.out.splitted_files, REFGEN_PARSEBIO.out)
    emit:
        DEMULTIPLEX.out.splitted_files
}
