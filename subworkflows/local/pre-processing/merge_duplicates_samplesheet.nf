//
// Subworkflow that merges FASTQ files of duplicate entries in a samplesheet
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MERGE_FASTQS } from '../../../modules/local/custom/manipulate/merge_files/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow merge_fastqs_subworkflow {

    take:
        ch_preproc_fastqs

        // Group by sample‑ID so we have all lanes / replicates for the same biological sample in one place
        ch_grouped = ch_preproc_fastqs                       \
            .map   { meta, fq1, fq2 -> tuple(meta.id, meta, fq1, fq2) } \
            .groupTuple()                                    \
            .map   { sid, rows ->
                        def meta        = rows[0][1]     // keep first meta block
                        def cDNA_list   = rows.collect { it[2] }.findAll{ it }  // drop nulls if any
                        def BCUMI_list  = rows.collect { it[3] }.findAll{ it }
                        tuple(meta, cDNA_list, BCUMI_list)
                    }

    main:
        // Concatenate all FASTQs belonging to the sample
        MERGE_FASTQS(ch_grouped)

    emit:
        data_output = MERGE_FASTQS.out
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
