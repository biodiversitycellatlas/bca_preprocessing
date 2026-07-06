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
include { SECONDDERIV_CELLCALLING                           } from '../../../modules/local/custom/dashboard/2nd_deriv/cellcalling/main'
include { FILTER_MATRICES                                   } from '../../../modules/local/custom/dashboard/2nd_deriv/filter_mtx/main'


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
        ref_fasta
        remap_geneext

    main:
        // Initialize reporting channels
        def ch_starsolo_bam     = Channel.empty()
        def ch_filtered_mtx     = Channel.empty()

        // Check if star index is provided, if not create it
        def star_index_ch
        if (params.star_index && file(params.star_index).exists() && remap_geneext == 'false') {
            star_index_ch = Channel.value(file(params.star_index))
        } else {
            STARSOLO_INDEX(data_output, ref_gtf, ref_fasta)
            star_index_ch = STARSOLO_INDEX.out.index.first()
        }

        // Confirm bc_whitelist is a safe value channel
        def bc_whitelist_ch = bc_whitelist instanceof List
            ? Channel.value(bc_whitelist)
            : bc_whitelist.ifEmpty([]).first()

        // Run STARsolo alignment
        STARSOLO_ALIGN(data_output, bc_whitelist, star_index_ch)

        // Filters the raw matrices based on the cutoff determined by the second derivative method, if selected
        if (params.cellfilter_method == "second_derivative") {
            SECONDDERIV_CELLCALLING(STARSOLO_ALIGN.out.umi_per_cell)
            FILTER_MATRICES(STARSOLO_ALIGN.out.genefull50_raw_dir, SECONDDERIV_CELLCALLING.out.cutoff)
            ch_filtered_mtx = FILTER_MATRICES.out.filtered_matrix
        } else {
            ch_filtered_mtx = STARSOLO_ALIGN.out.genefull50_filtered_dir
        }

        if (params.star_generateBAM) {
            ch_starsolo_bam = STARSOLO_ALIGN.out.bam_file
        }


    emit:
        mapping_files                   = STARSOLO_ALIGN.out.starsolo_files
        starsolo_bam                    = ch_starsolo_bam
        star_solodir                    = STARSOLO_ALIGN.out.star_solodir
        starsolo_genefull50_raw         = STARSOLO_ALIGN.out.genefull50_raw_dir
        starsolo_genefull50_filtered    = ch_filtered_mtx
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
