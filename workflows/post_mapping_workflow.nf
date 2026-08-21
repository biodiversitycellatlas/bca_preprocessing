//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { restage_mapping_workflow                                          } from '../subworkflows/local/post-processing/restage_mapping'
include { cellcalling_starsolo_workflow                                     } from '../subworkflows/local/post-processing/cellcalling'
include { cellcalling_alevin_workflow                                       } from '../subworkflows/local/post-processing/cellcalling'
include { bam_inspection_workflow                                           } from '../subworkflows/local/post-processing/bam_inspection'
include { bam_inspection_workflow as bam_inspection_geneext_workflow        } from '../subworkflows/local/post-processing/bam_inspection'

include { MERGE_REF_GTF                                                     } from '../modules/local/custom/manipulate/merge_ref_gtf/main'
include { MERGE_REF_GTF as MERGE_REF_GTF_GENEEXT                            } from '../modules/local/custom/manipulate/merge_ref_gtf/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RE-RUN EVERYTHING AFTER MAPPING
        Picks a previous run up at the point mapping ended: its published STARsolo and
        alevin-fry results are read back in, the cells are called again from the current
        samplesheet, and every step that depends on the cell set is redone on top --
        matrix filtering, the saturation curve, per-cell metrics, ambient-RNA removal,
        doublet detection and the report.

        Mapping itself, including the GeneExt re-mapping, is not repeated: when
        'perform_geneext' is set, the GeneExt-remapped analytical runs of the previous
        run are picked up alongside the standard ones and carried through the same steps.

        Emits the same channels as QC_mapping_workflow, so the filtering and reporting
        workflows are driven identically in both modes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow post_mapping_workflow {
    take:
        ch_samplesheet

    main:
        // Initialize reporting channels
        def ch_sat_imgs      = Channel.empty()
        def ch_sat_res_imgs  = Channel.empty()
        def ch_sat_logs      = Channel.empty()
        def ch_featurecounts = Channel.empty()
        def ch_pavian_sankey = Channel.empty()
        def ch_geneext_report = Channel.empty()
        def ch_geneext_log    = Channel.empty()

        def prev_dir = (params.previous_outdir ?: params.outdir).toString()

        // Read the previous run's mapping results back into channels
        restage_mapping_workflow(ch_samplesheet)

        // Conditionally bypass MERGE_REF_GTF when no additional features are provided
        def ref_gtf_ch
        if (params.ref_gtf_addfeature) {
            MERGE_REF_GTF(params.ref_gtf, Channel.fromPath(params.ref_gtf_addfeature))
            ref_gtf_ch = MERGE_REF_GTF.out.gtf
        } else {
            ref_gtf_ch = Channel.value(file(params.ref_gtf))
        }

        // The GeneExt-remapped runs were annotated against the extended GTF, so their
        // read-level statistics have to be recomputed against that same annotation
        def ref_gtf_geneext_ch = ref_gtf_ch
        if (params.perform_geneext) {
            def geneext_gtf = file(params.geneext_gtf ?: "${prev_dir}/gene_ext/geneext.gtf")
            if (!geneext_gtf.exists()) {
                error(
                    "'perform_geneext' is set, but no GeneExt annotation was found at ${geneext_gtf}.\n" +
                    "Point 'geneext_gtf' at the extended GTF of the run being picked up,\n" +
                    "or unset 'perform_geneext' to process the standard mapping runs only."
                )
            }
            if (params.ref_gtf_addfeature) {
                MERGE_REF_GTF_GENEEXT(Channel.value(geneext_gtf), Channel.fromPath(params.ref_gtf_addfeature))
                ref_gtf_geneext_ch = MERGE_REF_GTF_GENEEXT.out.gtf
            } else {
                ref_gtf_geneext_ch = Channel.value(geneext_gtf)
            }

            // Check if GeneExt produced a report and/or log, and if so, publish them to the dashboard
            def geneext_report = file("${geneext_gtf}.Report.html")
            def geneext_log    = file("${geneext_gtf}.GeneExt.log")
            if (geneext_report.exists()) {
                ch_geneext_report = Channel.value(geneext_report)
            }
            if (geneext_log.exists()) {
                ch_geneext_log = Channel.value(geneext_log)
            }
        }

        // Re-call cells. One invocation covers the standard and the GeneExt-remapped runs,
        // since the cutoff is derived per sample from that run's own UMI curve.
        cellcalling_starsolo_workflow(
            restage_mapping_workflow.out.starsolo_genefull50_raw,
            restage_mapping_workflow.out.star_umipercell,
            restage_mapping_workflow.out.star_cellreads,
            restage_mapping_workflow.out.starsolo_genefull50_filtered
        )

        cellcalling_alevin_workflow(restage_mapping_workflow.out.af_mtx)

        // Both mappers emit the same second-derivative artefacts, so the reporting channels carry them together
        def ch_secondderiv_knee = cellcalling_starsolo_workflow.out.secondderiv_knee
            .mix(cellcalling_alevin_workflow.out.secondderiv_knee)
        def ch_secondderiv_stats = cellcalling_starsolo_workflow.out.secondderiv_stats
            .mix(cellcalling_alevin_workflow.out.secondderiv_stats)
        def ch_secondderiv_cutoff = cellcalling_starsolo_workflow.out.secondderiv_cutoff
            .mix(cellcalling_alevin_workflow.out.secondderiv_cutoff)

        // The saturation curve is fitted per cell and featureCounts is annotation-specific,
        // so the two annotations are inspected separately, as during mapping
        if (params.star_generateBAM) {

            def ch_bam       = restage_mapping_workflow.out.starsolo_bam
            def ch_summary   = restage_mapping_workflow.out.star_summaries
            def ch_final_log = restage_mapping_workflow.out.star_final_log
            def ch_sd_stats  = cellcalling_starsolo_workflow.out.secondderiv_stats

            bam_inspection_workflow(
                ch_bam.filter       { meta, _bam  -> !meta.geneext },
                ref_gtf_ch,
                ch_summary.filter   { meta, _csv  -> !meta.geneext },
                ch_final_log.filter { meta, _log  -> !meta.geneext },
                ch_sd_stats.filter  { meta, _json -> !meta.geneext }
            )

            ch_sat_imgs      = bam_inspection_workflow.out.saturation_imgs
            ch_sat_res_imgs  = bam_inspection_workflow.out.saturation_residual_imgs
            ch_sat_logs      = bam_inspection_workflow.out.saturation_logs
            ch_featurecounts = bam_inspection_workflow.out.featurecount_txt
            ch_pavian_sankey = bam_inspection_workflow.out.pavian_sankey

            if (params.perform_geneext) {
                bam_inspection_geneext_workflow(
                    ch_bam.filter       { meta, _bam  -> meta.geneext },
                    ref_gtf_geneext_ch,
                    ch_summary.filter   { meta, _csv  -> meta.geneext },
                    ch_final_log.filter { meta, _log  -> meta.geneext },
                    ch_sd_stats.filter  { meta, _json -> meta.geneext }
                )

                ch_sat_imgs      = ch_sat_imgs.mix(bam_inspection_geneext_workflow.out.saturation_imgs)
                ch_sat_res_imgs  = ch_sat_res_imgs.mix(bam_inspection_geneext_workflow.out.saturation_residual_imgs)
                ch_sat_logs      = ch_sat_logs.mix(bam_inspection_geneext_workflow.out.saturation_logs)
                ch_featurecounts = ch_featurecounts.mix(bam_inspection_geneext_workflow.out.featurecount_txt)
                ch_pavian_sankey = ch_pavian_sankey.mix(bam_inspection_geneext_workflow.out.pavian_sankey)
            }
        }

    emit:
        mapped_samplesheet           = restage_mapping_workflow.out.mapped_metas
        ref_gtf                      = ref_gtf_ch
        mapping_files                = Channel.empty()
        starsolo_bam                 = restage_mapping_workflow.out.starsolo_bam
        star_solodir                 = restage_mapping_workflow.out.star_solodir
        starsolo_genefull50_raw      = restage_mapping_workflow.out.starsolo_genefull50_raw
        starsolo_genefull50_filtered = cellcalling_starsolo_workflow.out.filtered_matrix
        secondderiv_knee             = ch_secondderiv_knee
        secondderiv_stats            = ch_secondderiv_stats
        secondderiv_cutoff           = ch_secondderiv_cutoff
        saturation_imgs              = ch_sat_imgs
        saturation_residual_imgs     = ch_sat_res_imgs
        saturation_logs              = ch_sat_logs
        star_umipercell              = restage_mapping_workflow.out.star_umipercell
        star_log                     = restage_mapping_workflow.out.star_log
        star_final_log               = restage_mapping_workflow.out.star_final_log
        star_summaries               = restage_mapping_workflow.out.star_summaries
        star_cellreads               = restage_mapping_workflow.out.star_cellreads
        af_meta_info                 = restage_mapping_workflow.out.af_meta_info
        af_quant_json                = restage_mapping_workflow.out.af_quant_json
        af_cell_meta                 = restage_mapping_workflow.out.af_cell_meta
        af_mtx                       = restage_mapping_workflow.out.af_mtx
        af_filtered_mtx              = cellcalling_alevin_workflow.out.filtered_matrix
        af_umipercell                = cellcalling_alevin_workflow.out.umi_per_cell
        featurecount_txt             = ch_featurecounts
        pavian_sankey                = ch_pavian_sankey
        geneext_report               = ch_geneext_report
        geneext_log                  = ch_geneext_log
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
