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
include { cellcalling_workflow                              } from '../post-processing/cellcalling'


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
        bam_only            // true for a pass mapped only for GeneExt: alignment and BAM, nothing else
        prebuilt_star_index // a STARsolo index another pass already built, or null to build one here

    main:
        // Initialize reporting channels
        def ch_starsolo_bam         = channel.empty()
        def ch_filtered_mtx         = channel.empty()
        def ch_velocyto_filtered    = channel.empty()
        def ch_secondderiv_knee     = channel.empty()
        def ch_secondderiv_stats    = channel.empty()
        def ch_secondderiv_cutoff   = channel.empty()
        def ch_versions             = channel.empty()

        // Check if star index is provided, if not create it
        def star_index_ch
        if (prebuilt_star_index != null) {
            star_index_ch = prebuilt_star_index
        } else if (params.star_index && remap_geneext == 'false') {
            // A given index that is not there is a mistake, not a request to build one: the
            // build takes hours and is published under outdir, not at 'star_index'
            if (!file(params.star_index).exists()) {
                error("'star_index' = '${params.star_index}' does not exist (relative paths are read from the launch directory, ${workflow.launchDir}).\nUnset it to build an index in this run.")
            }
            star_index_ch = channel.value(file(params.star_index))
        } else {
            // Build a single index for the whole run; the samples' cDNA reads are only
            // gathered so the index can size --sjdbOverhang on the longest read
            def all_cdna_ch = data_output
                .map { meta, fastq_cDNA, fastq_BC_UMI, fastq_indices, input_file -> fastq_cDNA }
                .collect()
            def index_suffix = remap_geneext == 'true' ? '_geneext' : ''
            STARSOLO_INDEX(all_cdna_ch, ref_gtf, ref_fasta, index_suffix)
            star_index_ch = STARSOLO_INDEX.out.index.first()
            ch_versions = ch_versions.mix(STARSOLO_INDEX.out.versions)
        }

        // The whitelist reaches every sample: a value channel holding the list of files
        def bc_whitelist_ch = bc_whitelist instanceof List ? channel.value(bc_whitelist) : bc_whitelist

        // Protocol 'seqspec' reads its barcode layout from the seqspec file, staged into the task
        def seqspec = params.protocol == 'seqspec' ? file(params.seqspec_file, checkIfExists: true) : []

        // Run STARsolo alignment
        STARSOLO_ALIGN(data_output, bc_whitelist_ch, bam_only, star_index_ch, seqspec)
        ch_versions = ch_versions.mix(STARSOLO_ALIGN.out.versions)

        // Re-call cells on a UMI cutoff, or keep STARsolo's own filtered matrix
        if (!bam_only) {
            cellcalling_workflow(
                'starsolo',
                STARSOLO_ALIGN.out.genefull50_raw_dir,
                STARSOLO_ALIGN.out.umi_per_cell,
                STARSOLO_ALIGN.out.cellreads_stats,
                STARSOLO_ALIGN.out.genefull50_filtered_dir,
                STARSOLO_ALIGN.out.velocyto_raw_dir
            )

            ch_filtered_mtx       = cellcalling_workflow.out.filtered_matrix
            ch_velocyto_filtered  = cellcalling_workflow.out.velocyto_filtered
            ch_secondderiv_knee   = cellcalling_workflow.out.secondderiv_knee
            ch_secondderiv_stats  = cellcalling_workflow.out.secondderiv_stats
            ch_secondderiv_cutoff = cellcalling_workflow.out.secondderiv_cutoff
            ch_versions           = ch_versions.mix(cellcalling_workflow.out.versions)
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
        star_index                      = star_index_ch
        versions                        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
