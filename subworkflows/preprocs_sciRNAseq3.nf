//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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
workflow sciRNAseq3_workflow {
    take:
        ch_samplesheet
    main:
        // Demultiplex each splitted R1 and R2 file to generate sample-specific fastq.gz files
        SCIROCKET_DEMUX(ch_samplesheet)

        // Trimming adapters and low-quality reads
        FASTP(SCIROCKET_DEMUX.out.samples)

    emit:
        data_output     = FASTP.out
        bc_whitelist    = SCIROCKET_DEMUX.out.bc_whitelist
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/