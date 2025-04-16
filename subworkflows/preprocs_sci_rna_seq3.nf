//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_DATA         } from '../modules/download'
include { FASTQ_SPLITTER        } from '../modules/fastqsplitter'
include { SCIROCKET_DEMUX       } from '../modules/scirocket_demux'
include { SCIROCKET_SEQS_GATHER } from '../modules/scirocket_seqs_gather'
include { SCIROCKET_SPL_GATHER  } from '../modules/scirocket_spl_gather'
include { FASTP                 } from '../modules/fastp'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR SCI-RNA-SEQ3 DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow sci_rna_seq3_workflow {
    take:
        sample_ids
    main:
        // Import the fastq files into the nf workdir using sym links to the original files
        DOWNLOAD_DATA(sample_ids)

        // Splitting fastq files into smaller chunks
        FASTQ_SPLITTER(sample_ids, DOWNLOAD_DATA.out)

        // Flatten the output of FASTQ_SPLITTER to create a list of tuples
        FASTQ_SPLITTER.out.split_files_ch
            .flatMap { sample_id, r1_files, r2_files ->
                def paired = []
                // Assume that both lists have the same number of files equal to params.scatter
                for( int i = 0; i < r1_files.size(); i++ ) {
                    paired << [ sample_id, i+1, r1_files[i], r2_files[i] ]
                }
                return paired
            }
            .set { demux_input_ch }

        // Demultiplex each splitted R1 and R2 file to generate sample-specific fastq.gz files
        SCIROCKET_DEMUX(demux_input_ch)

        // Collecting the demultiplexed fastq files
        demux_scatter_ch = SCIROCKET_DEMUX.out.demux_scatter_ch.groupTuple()

        // Combining the discarded and whitelist files
        SCIROCKET_SEQS_GATHER(demux_scatter_ch)

        // Combining the sample-specific fastq.fz files
        SCIROCKET_SPL_GATHER(SCIROCKET_SEQS_GATHER.out.gathered_seq_ch)

        // Trimming adapters and low-quality reads
        FASTP(SCIROCKET_SPL_GATHER.out.final_samples_ch)

    emit:
        FASTP.out
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/