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
include { FASTP                 } from '../../../modules/local/tools/fastp/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR SCI-RNA-SEQ3 DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow sciRNAseq3_nogather_workflow {
    take:
        ch_samplesheet

    main:

        // Run demux per‐sample and capture its output channels
        SCIROCKET_DEMUX(ch_samplesheet)

        SCIROCKET_SEQS_GATHER(
            SCIROCKET_DEMUX.out.samples_discarded_R1,
            SCIROCKET_DEMUX.out.samples_discarded_R2,
            SCIROCKET_DEMUX.out.bc_whitelists_p5,
            SCIROCKET_DEMUX.out.bc_whitelists_p7,
            SCIROCKET_DEMUX.out.bc_whitelists_ligation,
            SCIROCKET_DEMUX.out.bc_whitelists_rt
        )

        // Key the samplesheet by its meta.id
        def ch_samplesheet_keyed = ch_samplesheet.map { rec ->
            tuple(rec.meta.id as String, rec)                 // -> (id, rec)
        }

        // Key DEMUX outputs by extracting ID from filenames
        def keyById = { Path p ->
            def base = p.getBaseName()                        // e.g. SAMPLEID_R1.fastq
            def id   = base.replaceFirst(/_R[12].*$/, '')     // -> SAMPLEID
            tuple(id as String, p)                            // -> (id, path)
        }
        def ch_r1_keyed = SCIROCKET_DEMUX.out.samples_R1.map(keyById)   // (id, r1)
        def ch_r2_keyed = SCIROCKET_DEMUX.out.samples_R2.map(keyById)   // (id, r2)

        // Join on id and re-wrap
        ch_wrapped = ch_samplesheet_keyed
            .join(ch_r1_keyed)                                // (id, rec, r1)
            .join(ch_r2_keyed)                                // (id, rec, r1, r2)
            .map { id, rec, r1, r2 ->
                [
                    meta         : rec.meta,                  // keep metadata from samplesheet
                    fastq_cDNA   : r2,                        // R2 = cDNA
                    fastq_BC_UMI : r1,                        // R1 = BC+UMI
                    input_file   : rec.input_file             // preserve original input_file
                ]
            }

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
