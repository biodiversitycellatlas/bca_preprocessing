//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { mapping_starsolo_workflow                                         } from '../subworkflows/local/mapping/mapping_starsolo'
include { mapping_starsolo_workflow as mapping_starsolo_geneext_workflow    } from '../subworkflows/local/mapping/mapping_starsolo'
include { mapping_alevin_workflow                                           } from '../subworkflows/local/mapping/mapping_alevin'
include { bam_inspection_workflow                                           } from '../subworkflows/local/post-processing/bam_inspection'
include { bam_inspection_workflow as bam_inspection_geneext_workflow        } from '../subworkflows/local/post-processing/bam_inspection'
include { geneext_workflow                                                  } from '../subworkflows/local/mapping/geneext'

include { FASTQC                    } from '../modules/local/tools/fastqc/main'
include { SUBSAMPLE_FASTQS          } from '../modules/local/tools/seqtk/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow QC_mapping_workflow {
    take:
        data_output
        bc_whitelist

    main:
        // Initialize reporting channels
        def ch_mapping_files = Channel.empty()
        def ch_starsolo_bam  = Channel.empty()
        def ch_star_solodir  = Channel.empty()
        def ch_starsolo_genefull50_raw  = Channel.empty()
        def ch_starsolo_genefull50_filtered  = Channel.empty()
        def ch_sat_imgs         = Channel.empty()
        def ch_sat_res_imgs     = Channel.empty()
        def ch_sat_logs         = Channel.empty()
        def ch_star_umi         = Channel.empty()
        def ch_star_log         = Channel.empty()
        def ch_star_final_log   = Channel.empty()
        def ch_star_summaries   = Channel.empty()
        def ch_star_cellreads  = Channel.empty()
        def ch_alevin_meta_info = Channel.empty()
        def ch_alevin_quant_json = Channel.empty()
        def ch_alevin_cell_meta = Channel.empty()
        def ch_alevin_mtx       = Channel.empty()
        def ch_featurecounts    = Channel.empty()
        def ch_pavian_sankey    = Channel.empty()

        // Quality Control
        if (params.protocol != "scalebio") {
            FASTQC(data_output)
        }

        // Add suffix to metadata for seperation of results
        def apply_suffix = { ch, suffix ->
            ch.map { meta, fq1, fq2, fqidx, ss ->
                def new_meta = meta.clone()
                new_meta.base_id = meta.id          // Preserve 'sample1'
                new_meta.id = meta.id + suffix      // Creates 'sample1_starsolo'
                [new_meta, fq1, fq2, fqidx, ss]
            }
        }

        // Mapping: starsolo, alevin, alevin_starsolo (both), or alevin_subsampled_starsolo
        if (params.mapping_software == "starsolo" || params.mapping_software == "both" || params.mapping_software == "alevin_subsampled_starsolo" || params.mapping_software == "alevin_starsolo") {

            // If 'alevin_subsampled_starsolo' is selected, run STARsolo on a subsampled dataset
            if (params.mapping_software == "alevin_subsampled_starsolo") {
                def ch_star = apply_suffix(data_output, "_subsampled_starsolo")

                SUBSAMPLE_FASTQS(ch_star)
                mapping_starsolo_workflow(SUBSAMPLE_FASTQS.out, bc_whitelist, params.ref_gtf, 'false')

            // Run STARsolo on the full dataset
            } else {
                def ch_star = apply_suffix(data_output, "_starsolo")
                mapping_starsolo_workflow(ch_star, bc_whitelist, params.ref_gtf, 'false')
            }

            // Assign starsolo standard mapping outputs
            ch_mapping_files         =  mapping_starsolo_workflow.out.mapping_files
            ch_starsolo_bam          =  mapping_starsolo_workflow.out.starsolo_bam
            ch_star_solodir          =  mapping_starsolo_workflow.out.star_solodir
            ch_starsolo_genefull50_raw  = mapping_starsolo_workflow.out.starsolo_genefull50_raw
            ch_starsolo_genefull50_filtered = mapping_starsolo_workflow.out.starsolo_genefull50_filtered
            ch_star_umi              =  mapping_starsolo_workflow.out.star_umipercell
            ch_star_log              =  mapping_starsolo_workflow.out.star_log
            ch_star_final_log        =  mapping_starsolo_workflow.out.star_final_log
            ch_star_summaries        =  mapping_starsolo_workflow.out.star_summaries
            ch_star_cellreads        =  mapping_starsolo_workflow.out.star_cellreads


            // Run BAM inspection workflow on standard STARsolo output
            if (params.star_generateBAM) {
                bam_inspection_workflow(ch_starsolo_bam, ch_star_summaries, ch_star_final_log)

                ch_sat_imgs              =  bam_inspection_workflow.out.saturation_imgs
                ch_sat_res_imgs          =  bam_inspection_workflow.out.saturation_residual_imgs
                ch_sat_logs              =  bam_inspection_workflow.out.saturation_logs
                ch_featurecounts         =  bam_inspection_workflow.out.featurecount_txt
                ch_pavian_sankey         =  bam_inspection_workflow.out.pavian_sankey
            }

            // Optionally run geneext and rerun mapping steps
            if (params.perform_geneext || params.run_method == "geneext_only") {
                geneext_workflow(mapping_starsolo_workflow.out.starsolo_bam)

                if (params.perform_geneext) {
                    def ch_geneext = apply_suffix(data_output, "_geneext_starsolo")   // TODO: add _geneext_subsampled_starsolo option
                    mapping_starsolo_geneext_workflow(ch_geneext, bc_whitelist, geneext_workflow.out.ref_gtf, 'true')

                    // Mix the geneext starsolo outputs with standard run
                    ch_mapping_files                = ch_mapping_files.mix(mapping_starsolo_geneext_workflow.out.mapping_files)
                    ch_starsolo_bam                 = ch_starsolo_bam.mix(mapping_starsolo_geneext_workflow.out.starsolo_bam)
                    ch_star_solodir                 = ch_star_solodir.mix(mapping_starsolo_geneext_workflow.out.star_solodir)
                    ch_starsolo_genefull50_raw      = ch_starsolo_genefull50_raw.mix(mapping_starsolo_geneext_workflow.out.starsolo_genefull50_raw)
                    ch_starsolo_genefull50_filtered = ch_starsolo_genefull50_filtered.mix(mapping_starsolo_geneext_workflow.out.starsolo_genefull50_filtered)
                    ch_star_umi                     = ch_star_umi.mix(mapping_starsolo_geneext_workflow.out.star_umipercell)
                    ch_star_log                     = ch_star_log.mix(mapping_starsolo_geneext_workflow.out.star_log)
                    ch_star_final_log               = ch_star_final_log.mix(mapping_starsolo_geneext_workflow.out.star_final_log)
                    ch_star_summaries               = ch_star_summaries.mix(mapping_starsolo_geneext_workflow.out.star_summaries)
                    ch_star_cellreads               = ch_star_cellreads.mix(mapping_starsolo_geneext_workflow.out.star_cellreads)

                    // Run BAM inspection workflow on geneext STARsolo output
                    if (params.star_generateBAM) {
                        bam_inspection_geneext_workflow(mapping_starsolo_geneext_workflow.out.starsolo_bam,
                                                    mapping_starsolo_geneext_workflow.out.star_summaries,
                                                    mapping_starsolo_geneext_workflow.out.star_final_log)

                        ch_sat_imgs                     = ch_sat_imgs.mix(bam_inspection_geneext_workflow.out.saturation_imgs)
                        ch_sat_res_imgs                 = ch_sat_res_imgs.mix(bam_inspection_geneext_workflow.out.saturation_residual_imgs)
                        ch_sat_logs                     = ch_sat_logs.mix(bam_inspection_geneext_workflow.out.saturation_logs)
                        ch_featurecounts                = ch_featurecounts.mix(bam_inspection_geneext_workflow.out.featurecount_txt)
                        ch_pavian_sankey                = ch_pavian_sankey.mix(bam_inspection_geneext_workflow.out.pavian_sankey)
                    }
                }
            }
        }

        if (params.mapping_software == "alevin" || params.mapping_software == "both" || params.mapping_software == "alevin_subsampled_starsolo" || params.mapping_software == "alevin_starsolo") {
            def ch_alevin = apply_suffix(data_output, "_alevinfry")
            mapping_alevin_workflow(ch_alevin, bc_whitelist)

            ch_mapping_files = ch_mapping_files.mix(mapping_alevin_workflow.out.mapping_files)
            ch_alevin_meta_info = mapping_alevin_workflow.out.af_meta_info
            ch_alevin_quant_json = mapping_alevin_workflow.out.af_quant_json
            ch_alevin_cell_meta = mapping_alevin_workflow.out.af_cell_meta
            ch_alevin_mtx = mapping_alevin_workflow.out.af_mtx

        }

    emit:
        mapping_files            = ch_mapping_files
        starsolo_bam             = ch_starsolo_bam
        star_solodir             = ch_star_solodir
        starsolo_genefull50_raw  = ch_starsolo_genefull50_raw
        starsolo_genefull50_filtered = ch_starsolo_genefull50_filtered
        saturation_imgs          = ch_sat_imgs
        saturation_residual_imgs = ch_sat_res_imgs
        saturation_logs          = ch_sat_logs
        star_umipercell          = ch_star_umi
        star_log                 = ch_star_log
        star_final_log           = ch_star_final_log
        star_summaries           = ch_star_summaries
        star_cellreads           = ch_star_cellreads
        af_meta_info             = ch_alevin_meta_info
        af_quant_json            = ch_alevin_quant_json
        af_cell_meta             = ch_alevin_cell_meta
        af_mtx                   = ch_alevin_mtx
        featurecount_txt         = ch_featurecounts
        pavian_sankey            = ch_pavian_sankey
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
