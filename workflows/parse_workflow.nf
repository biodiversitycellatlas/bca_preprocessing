include { DOWNLOAD_DATA } from '../modules/download'
include { DEMUX_SPIPE } from '../modules/demux_spipe'
include { DEMUX_UMITOOLS } from '../modules/demux_umitools_parsebio'
include { REFGEN_PARSEBIO } from '../modules/refgen_parsebio'
include { MAPPING_PARSEBIO } from '../modules/mapping_parsebio'


/*
 * Parse Biosciences workflow
 *
 * Based on the sample wells, the fastq sequences must be split 
 * based on the first barcode round. Performs both custom- and 
 * commercial pre-processing (using the Parse Biosciences pipeline v1.3.1) 
 * enabling comparison between the methods and validation of steps.   
 */    

workflow parse_workflow {
    take:
        sample_ids
    main:
        // Import the fastq files into the nf workdir using sym links to the original files
        DOWNLOAD_DATA(sample_ids)
        
        // For each sample_id, all groups are added to the Channel
        groups = Channel.fromList(params.groups)
        comb_data = DOWNLOAD_DATA.out.fastq_files.combine(groups)

        // Demultiplex the fastq files based on the sample wells
        DEMUX_SPIPE(comb_data)
        DEMUX_UMITOOLS_PARSEBIO(comb_data)
        
        // Run the Parse Biosciences pipeline
        REFGEN_PARSEBIO()
        MAPPING_PARSEBIO(DEMUX_SPIPE.out.splitted_files, REFGEN_PARSEBIO.out)
    
    emit:
        // Result: demultiplexed fastq files from split-pipe
        DEMUX_SPIPE.out.splitted_files
}
