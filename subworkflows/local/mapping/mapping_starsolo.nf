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
include { cellcalling_starsolo_workflow                     } from '../post-processing/cellcalling'


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
        def ch_starsolo_bam         = Channel.empty()
        def ch_filtered_mtx         = Channel.empty()
        def ch_velocyto_filtered    = Channel.empty()
        def ch_secondderiv_knee     = Channel.empty()
        def ch_secondderiv_stats    = Channel.empty()
        def ch_secondderiv_cutoff   = Channel.empty()

        // Check if the user has requested a valid geneext-only run
        def bam_only = params.geneext_bam_only && params.run_method == "geneext_only" && remap_geneext == 'false'

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
        STARSOLO_ALIGN(data_output, bc_whitelist, bam_only, star_index_ch)

        // Re-call cells on a UMI cutoff, or keep STARsolo's own filtered matrix
        if (!bam_only) {
            cellcalling_starsolo_workflow(
                STARSOLO_ALIGN.out.genefull50_raw_dir,
                STARSOLO_ALIGN.out.umi_per_cell,
                STARSOLO_ALIGN.out.cellreads_stats,
                STARSOLO_ALIGN.out.genefull50_filtered_dir,
                STARSOLO_ALIGN.out.velocyto_raw_dir
            )

            ch_filtered_mtx       = cellcalling_starsolo_workflow.out.filtered_matrix
            ch_velocyto_filtered  = cellcalling_starsolo_workflow.out.velocyto_filtered
            ch_secondderiv_knee   = cellcalling_starsolo_workflow.out.secondderiv_knee
            ch_secondderiv_stats  = cellcalling_starsolo_workflow.out.secondderiv_stats
            ch_secondderiv_cutoff = cellcalling_starsolo_workflow.out.secondderiv_cutoff
        }

        // 'bam_only' forces the BAM on regardless of star_generateBAM: it is the only
        // output the mode produces, and GeneExt would otherwise be handed nothing
        if (params.star_generateBAM || bam_only) {
            ch_starsolo_bam = STARSOLO_ALIGN.out.bam_file
        }


    emit:
        mapping_files                   = STARSOLO_ALIGN.out.starsolo_files
        starsolo_bam                    = ch_starsolo_bam
        star_solodir                    = STARSOLO_ALIGN.out.star_solodir
        starsolo_genefull50_raw         = STARSOLO_ALIGN.out.genefull50_raw_dir
        starsolo_genefull50_filtered    = ch_filtered_mtx
        starsolo_velocyto_raw           = STARSOLO_ALIGN.out.velocyto_raw_dir
        starsolo_velocyto_filtered      = ch_velocyto_filtered
        secondderiv_knee                = ch_secondderiv_knee
        secondderiv_stats               = ch_secondderiv_stats
        secondderiv_cutoff              = ch_secondderiv_cutoff
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
