//
// Subworkflow with functionality specific to the workflow 'mapping_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { STARSOLO_INDEX                                    } from '../../../modules/local/tools/star/starsolo_genome_generate/main'
include { STARSOLO_ALIGN                                    } from '../../../modules/local/tools/star/starsolo_align/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN STARSOLO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow mapping_starsolo_workflow {
    take:
        data_output
        bc_whitelist
        ref_gtf
        remap_geneext

    main:
        // Initialize reporting channels
        def ch_starsolo_bam  = Channel.empty()

        // Get the first cDNA fastq file from the samplesheet to use for STAR index generation
        def ch_first_cDNA = data_output
            .map { meta, fastq_cDNA, fastq_BC_UMI, empty_list, samplesheet -> fastq_cDNA }
            .first()

        // Check if star index is provided, if not create it
        def star_index_ch
        if (params.star_index && file(params.star_index).exists() && remap_geneext == 'false') {
            star_index_ch = Channel.value(file(params.star_index))
        } else {
            STARSOLO_INDEX(ref_gtf, ch_first_cDNA)
            star_index_ch = STARSOLO_INDEX.out
        }

        // Run STARsolo alignment
        STARSOLO_ALIGN(data_output, bc_whitelist, star_index_ch)

        if (params.star_generateBAM) {
            ch_starsolo_bam = STARSOLO_ALIGN.out.bam_file
        }


    emit:
        mapping_files                   = STARSOLO_ALIGN.out.starsolo_files
        starsolo_bam                    = ch_starsolo_bam
        star_solodir                    = STARSOLO_ALIGN.out.star_solodir
        starsolo_genefull50_raw         = STARSOLO_ALIGN.out.genefull50_raw_dir
        starsolo_genefull50_filtered    = STARSOLO_ALIGN.out.genefull50_filtered_dir
        star_umipercell                 = STARSOLO_ALIGN.out.umi_per_cell
        star_log                        = STARSOLO_ALIGN.out.log_file
        star_final_log                  = STARSOLO_ALIGN.out.log_final_file
        star_summaries                  = STARSOLO_ALIGN.out.summary_csv
        star_cellreads                  = STARSOLO_ALIGN.out.cellreads_stats
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
