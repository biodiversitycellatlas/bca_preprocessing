//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SCIROCKET_DEMUX       } from '../../../modules/local/tools/scirocket/scirocket_demux/main'
include { SCIROCKET_SEQS_GATHER } from '../../../modules/local/tools/scirocket/scirocket_seqs_gather/main'
include { SCIROCKET_SPL_GATHER  } from '../../../modules/local/tools/scirocket/scirocket_spl_gather/main'
include { FASTP                 } from '../../../modules/local/tools/fastp/main'

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
        all_samples_R1 = SCIROCKET_DEMUX.out.samples_R1.collect()
        all_samples_R2 = SCIROCKET_DEMUX.out.samples_R2.collect()
        all_discarded_R1  = SCIROCKET_DEMUX.out.samples_discarded_R1.collect()
        all_discarded_R2  = SCIROCKET_DEMUX.out.samples_discarded_R2.collect()
        all_whitelists_p5 = SCIROCKET_DEMUX.out.bc_whitelists_p5.collect()
        all_whitelists_p7 = SCIROCKET_DEMUX.out.bc_whitelists_p7.collect()
        all_whitelists_ligation = SCIROCKET_DEMUX.out.bc_whitelists_ligation.collect()
        all_whitelists_rt = SCIROCKET_DEMUX.out.bc_whitelists_rt.collect()

        // Runs each gather once after every demux has finished
        SCIROCKET_SEQS_GATHER(
            all_discarded_R1,
            all_discarded_R2,
            all_whitelists_p5,
            all_whitelists_p7,
            all_whitelists_ligation,
            all_whitelists_rt
        )
        SCIROCKET_SPL_GATHER(all_samples_R1, all_samples_R2)

        // Re-wrap the meta so FASTP sees meta.id
        wrapped = SCIROCKET_SPL_GATHER.out
            .map { sampleId, r1_file, r2_file ->
                def meta = [
                    id                : sampleId,
                    expected_cells    : 3000,
                ]

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
