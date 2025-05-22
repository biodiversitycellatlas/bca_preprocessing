//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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

        // Run demux per‐sample and capture its output channels
        SCIROCKET_DEMUX(ch_samplesheet)

        // Collect all of each channel into a single list
        all_discarded   = SCIROCKET_DEMUX.out.samples_discarded.collect()
        all_whitelists  = SCIROCKET_DEMUX.out.bc_whitelists.collect()
        // all_samples     = SCIROCKET_DEMUX.out.samples.collect()

        // Flatten the per‐sample lists into a single stream of FASTQ paths
        ch_all_demux = SCIROCKET_DEMUX.out.samples.flatten()

        // From samples channel, pull out R1 vs R2
        ch_all_R1 = ch_all_demux
                    .filter { it.name.endsWith('_R1.fastq.gz') }
                    .collect()
        ch_all_R2 = ch_all_demux
                    .filter { it.name.endsWith('_R2.fastq.gz') }
                    .collect()

        // Runs each gather once after every demux has finished
        SCIROCKET_SEQS_GATHER(all_discarded, all_whitelists)
        SCIROCKET_SPL_GATHER(ch_all_R1, ch_all_R2)

        // Re-wrap the meta so FASTP sees meta.id
        wrapped = SCIROCKET_SPL_GATHER.out
            .map { sampleId, r1_file, r2_file ->
                def meta         = [ id: sampleId ]

                // Return a named map
                [
                    meta         : meta,
                    fastq_cDNA   : r2_file,
                    fastq_BC_UMI : r1_file
                ]
            }
            .set { ch_wrapped }

        // Trimming adapters and low-quality reads
        FASTP(ch_wrapped)

    emit:
        data_output     = FASTP.out
        bc_whitelist    = SCIROCKET_SEQS_GATHER.out.bc_whitelist
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/