include { DOWNLOAD_DATA } from '../modules/download'
include { DEMUX_UMITOOLS_PARSEBIO } from '../modules/demux_umitools_parsebio'

include { PARSEBIO_PIPELINE_DEMUX } from '../modules/parsebio_pipeline_demux'
include { PARSEBIO_PIPELINE_MKREF } from '../modules/parsebio_pipeline_mkref'
include { PARSEBIO_PIPELINE } from '../modules/parsebio_pipeline'


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
        PARSEBIO_PIPELINE_DEMUX(comb_data)
        DEMUX_UMITOOLS_PARSEBIO(comb_data)
        
        // Run the Parse Biosciences pipeline
        PARSEBIO_PIPELINE_MKREF()
        PARSEBIO_PIPELINE(PARSEBIO_PIPELINE_DEMUX.out.splitted_files, PARSEBIO_PIPELINE_MKREF.out)
    
    emit:
        // Result: demultiplexed fastq files from split-pipe
        PARSEBIO_PIPELINE_DEMUX.out.splitted_files
}
