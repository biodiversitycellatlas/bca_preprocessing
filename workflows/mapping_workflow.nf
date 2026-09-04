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

include { MERGE_REF_FASTA                                                   } from '../modules/local/custom/manipulate/merge_ref_fasta/main'
include { MERGE_REF_FASTA as MERGE_REF_FASTA_GENEEXT                        } from '../modules/local/custom/manipulate/merge_ref_fasta/main'
include { MERGE_REF_GTF                                                     } from '../modules/local/custom/manipulate/merge_ref_gtf/main'
include { MERGE_REF_GTF as MERGE_REF_GTF_GENEEXT                            } from '../modules/local/custom/manipulate/merge_ref_gtf/main'
include { FASTQC                                                            } from '../modules/local/tools/fastqc/main'
include { SUBSAMPLE_FASTQS                                                  } from '../modules/local/tools/seqtk/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Name every row for the analytical run it belongs to, keeping the sample's own id in
// 'base_id'. Uses index access to avoid fixed-arity destructuring, since the rows carry
// different numbers of file elements.
def apply_suffix(ch, String suffix) {
    return ch.map { row ->
        def new_meta = row[0].clone()
        new_meta.base_id = row[0].id
        new_meta.id = row[0].id + suffix
        [new_meta] + row.drop(1)
    }
}

workflow QC_mapping_workflow {
    take:
        data_output
        bc_whitelist

    main:
        // Initialize channels
        def ch_samples = data_output
        def ch_mapped_ss                    = Channel.empty()
        def ch_mapping_files                = Channel.empty()
        def ch_starsolo_bam                 = Channel.empty()
        def ch_star_solodir                 = Channel.empty()
        def ch_starsolo_genefull50_raw      = Channel.empty()
        def ch_starsolo_genefull50_filtered = Channel.empty()
        def ch_starsolo_velocyto_raw        = Channel.empty()
        def ch_starsolo_velocyto_filtered   = Channel.empty()
        def ch_secondderiv_knee             = Channel.empty()
        def ch_secondderiv_stats            = Channel.empty()
        def ch_secondderiv_cutoff           = Channel.empty()
        def ch_sat_imgs                     = Channel.empty()
        def ch_sat_res_imgs                 = Channel.empty()
        def ch_sat_logs                     = Channel.empty()
        def ch_star_umi                     = Channel.empty()
        def ch_star_log                     = Channel.empty()
        def ch_star_final_log               = Channel.empty()
        def ch_star_summaries               = Channel.empty()
        def ch_star_cellreads               = Channel.empty()
        def ch_alevin_meta_info             = Channel.empty()
        def ch_alevin_quant_json            = Channel.empty()
        def ch_alevin_cell_meta             = Channel.empty()
        def ch_alevin_mtx                   = Channel.empty()
        def ch_alevin_filtered_mtx          = Channel.empty()
        def ch_alevin_umipercell            = Channel.empty()
        def ch_featurecounts                = Channel.empty()
        def ch_pavian_sankey                = Channel.empty()
        def ch_geneext_report               = Channel.empty()
        def ch_geneext_log                  = Channel.empty()

        // Switch to disable non-GeneExt outputs when the user has requested a geneext_only run
        def geneext_bam_only = params.geneext_bam_only && params.run_method == "geneext_only"

        // Conditionally bypass MERGE_REF_GTF/FASTA when no additional features are provided
        def ref_gtf_ch
        if (params.ref_gtf_addfeature) {
            MERGE_REF_GTF(params.ref_gtf, Channel.fromPath(params.ref_gtf_addfeature))
            ref_gtf_ch = MERGE_REF_GTF.out.gtf
        } else {
            ref_gtf_ch = Channel.value(file(params.ref_gtf))
        }

        def ref_fasta_ch
        if (params.ref_fasta_addfeature) {
            MERGE_REF_FASTA(params.ref_fasta, Channel.fromPath(params.ref_fasta_addfeature))
            ref_fasta_ch = MERGE_REF_FASTA.out.fasta
        } else {
            ref_fasta_ch = Channel.value(file(params.ref_fasta))
        }

        // Safe bc_whitelist: emit empty string when no whitelist is produced by preprocessing
        def bc_whitelist_safe = bc_whitelist.ifEmpty("")

        // Alevin takes the whitelist as a path input, which cannot be an empty string, so protocols running without a whitelist (e.g. MARS-seq) stage nothing instead
        def bc_whitelist_alevin = bc_whitelist_safe
            .map { wl -> wl?.toString()?.trim() ? wl : [] }

        // Quality Control
        FASTQC(ch_samples)

        // Mapping: starsolo, alevin, alevin_starsolo (both), or alevin_subsampled_starsolo
        if (params.mapping_software == "starsolo" || params.mapping_software == "both" || params.mapping_software == "alevin_subsampled_starsolo" || params.mapping_software == "alevin_starsolo") {

            // If 'alevin_subsampled_starsolo' is selected, run STARsolo on a subsampled dataset
            if (params.mapping_software == "alevin_subsampled_starsolo") {

                def ch_star           = apply_suffix(ch_samples, "_subsampled_starsolo")
                ch_mapped_ss = ch_mapped_ss.mix(ch_star)
                SUBSAMPLE_FASTQS(ch_star)

                mapping_starsolo_workflow(SUBSAMPLE_FASTQS.out.subsampled_files, bc_whitelist_safe, ref_gtf_ch, ref_fasta_ch, 'false')

            // Run STARsolo on the full dataset
            } else {

                def ch_star           = apply_suffix(ch_samples, "_starsolo")
                ch_mapped_ss = ch_mapped_ss.mix(ch_star)

                mapping_starsolo_workflow(ch_star, bc_whitelist_safe, ref_gtf_ch, ref_fasta_ch, 'false')
            }

            // Assign starsolo standard mapping outputs
            ch_mapping_files             =  mapping_starsolo_workflow.out.mapping_files
            ch_starsolo_bam              =  mapping_starsolo_workflow.out.starsolo_bam
            ch_star_solodir              =  mapping_starsolo_workflow.out.star_solodir
            ch_starsolo_genefull50_raw   =  mapping_starsolo_workflow.out.starsolo_genefull50_raw
            ch_starsolo_genefull50_filtered = mapping_starsolo_workflow.out.starsolo_genefull50_filtered
            ch_starsolo_velocyto_raw     =  mapping_starsolo_workflow.out.starsolo_velocyto_raw
            ch_starsolo_velocyto_filtered = mapping_starsolo_workflow.out.starsolo_velocyto_filtered
            ch_secondderiv_knee          =  mapping_starsolo_workflow.out.secondderiv_knee
            ch_secondderiv_stats         =  mapping_starsolo_workflow.out.secondderiv_stats
            ch_secondderiv_cutoff        =  mapping_starsolo_workflow.out.secondderiv_cutoff
            ch_star_umi                  =  mapping_starsolo_workflow.out.star_umipercell
            ch_star_log                  =  mapping_starsolo_workflow.out.star_log
            ch_star_final_log            =  mapping_starsolo_workflow.out.star_final_log
            ch_star_summaries            =  mapping_starsolo_workflow.out.star_summaries
            ch_star_cellreads            =  mapping_starsolo_workflow.out.star_cellreads

            // Run BAM inspection workflow on STARsolo output when run_method is not set to 'geneext_only'
            if (params.star_generateBAM && !geneext_bam_only) {
                bam_inspection_workflow(ch_starsolo_bam, ref_gtf_ch, ch_star_summaries, ch_star_final_log,
                                            mapping_starsolo_workflow.out.secondderiv_stats)

                ch_sat_imgs              =  bam_inspection_workflow.out.saturation_imgs
                ch_sat_res_imgs          =  bam_inspection_workflow.out.saturation_residual_imgs
                ch_sat_logs              =  bam_inspection_workflow.out.saturation_logs
                ch_featurecounts         =  bam_inspection_workflow.out.featurecount_txt
                ch_pavian_sankey         =  bam_inspection_workflow.out.pavian_sankey
            }

            // Optionally run geneext and rerun mapping steps
            if (params.perform_geneext || params.run_method == "geneext_only") {

                // Collect ALL BAMs from all samples before running geneext
                geneext_workflow(mapping_starsolo_workflow.out.starsolo_bam)

                // GeneExt's own run statistics, summarised in the dashboard
                ch_geneext_report = geneext_workflow.out.report
                ch_geneext_log    = geneext_workflow.out.geneext_log

                if (params.perform_geneext) {

                    // Same conditional bypass for geneext GTF/FASTA
                    def ref_gtf_geneext_ch
                    if (params.ref_gtf_addfeature) {
                        MERGE_REF_GTF_GENEEXT(geneext_workflow.out.ref_gtf, Channel.fromPath(params.ref_gtf_addfeature))
                        ref_gtf_geneext_ch = MERGE_REF_GTF_GENEEXT.out.gtf
                    } else {
                        // Geneext always extends from the geneext output, no bypass possible here
                        MERGE_REF_GTF_GENEEXT(geneext_workflow.out.ref_gtf, Channel.value([]))
                        ref_gtf_geneext_ch = MERGE_REF_GTF_GENEEXT.out.gtf
                    }

                    def ref_fasta_geneext_ch
                    if (params.ref_fasta_addfeature) {
                        MERGE_REF_FASTA_GENEEXT(params.ref_fasta, Channel.fromPath(params.ref_fasta_addfeature))
                        ref_fasta_geneext_ch = MERGE_REF_FASTA_GENEEXT.out.fasta
                    } else {
                        ref_fasta_geneext_ch = Channel.value(file(params.ref_fasta))
                    }

                    def ch_geneext           = apply_suffix(ch_samples, "_geneext_starsolo")
                    ch_mapped_ss = ch_mapped_ss.mix(ch_geneext)

                    mapping_starsolo_geneext_workflow(ch_geneext, bc_whitelist_safe, ref_gtf_geneext_ch, ref_fasta_geneext_ch, 'true')

                    // Mix the geneext starsolo outputs with standard run
                    ch_mapping_files                = ch_mapping_files.mix(mapping_starsolo_geneext_workflow.out.mapping_files)
                    ch_starsolo_bam                 = ch_starsolo_bam.mix(mapping_starsolo_geneext_workflow.out.starsolo_bam)
                    ch_star_solodir                 = ch_star_solodir.mix(mapping_starsolo_geneext_workflow.out.star_solodir)
                    ch_starsolo_genefull50_raw      = ch_starsolo_genefull50_raw.mix(mapping_starsolo_geneext_workflow.out.starsolo_genefull50_raw)
                    ch_starsolo_genefull50_filtered = ch_starsolo_genefull50_filtered.mix(mapping_starsolo_geneext_workflow.out.starsolo_genefull50_filtered)
                    ch_starsolo_velocyto_raw        = ch_starsolo_velocyto_raw.mix(mapping_starsolo_geneext_workflow.out.starsolo_velocyto_raw)
                    ch_starsolo_velocyto_filtered   = ch_starsolo_velocyto_filtered.mix(mapping_starsolo_geneext_workflow.out.starsolo_velocyto_filtered)
                    ch_secondderiv_knee             = ch_secondderiv_knee.mix(mapping_starsolo_geneext_workflow.out.secondderiv_knee)
                    ch_secondderiv_stats            = ch_secondderiv_stats.mix(mapping_starsolo_geneext_workflow.out.secondderiv_stats)
                    ch_secondderiv_cutoff           = ch_secondderiv_cutoff.mix(mapping_starsolo_geneext_workflow.out.secondderiv_cutoff)
                    ch_star_umi                     = ch_star_umi.mix(mapping_starsolo_geneext_workflow.out.star_umipercell)
                    ch_star_log                     = ch_star_log.mix(mapping_starsolo_geneext_workflow.out.star_log)
                    ch_star_final_log               = ch_star_final_log.mix(mapping_starsolo_geneext_workflow.out.star_final_log)
                    ch_star_summaries               = ch_star_summaries.mix(mapping_starsolo_geneext_workflow.out.star_summaries)
                    ch_star_cellreads               = ch_star_cellreads.mix(mapping_starsolo_geneext_workflow.out.star_cellreads)

                    // Run BAM inspection workflow on geneext STARsolo output
                    if (params.star_generateBAM) {
                        bam_inspection_geneext_workflow(mapping_starsolo_geneext_workflow.out.starsolo_bam,
                                                    ref_gtf_geneext_ch,
                                                    mapping_starsolo_geneext_workflow.out.star_summaries,
                                                    mapping_starsolo_geneext_workflow.out.star_final_log,
                                                    mapping_starsolo_geneext_workflow.out.secondderiv_stats)

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
            def ch_alevin = apply_suffix(ch_samples, "_alevinfry")
            ch_mapped_ss = ch_mapped_ss.mix(ch_alevin)
            mapping_alevin_workflow(ch_alevin, bc_whitelist_alevin)

            ch_mapping_files    = ch_mapping_files.mix(mapping_alevin_workflow.out.mapping_files)
            ch_alevin_meta_info = mapping_alevin_workflow.out.af_meta_info
            ch_alevin_quant_json = mapping_alevin_workflow.out.af_quant_json
            ch_alevin_cell_meta = mapping_alevin_workflow.out.af_cell_meta
            ch_alevin_mtx       = mapping_alevin_workflow.out.af_mtx
            ch_alevin_filtered_mtx = mapping_alevin_workflow.out.af_filtered_mtx
            ch_alevin_umipercell = mapping_alevin_workflow.out.af_umipercell

            // Both mappers emit the same second-derivative artefacts, so the reporting channels carry them together
            ch_secondderiv_knee   = ch_secondderiv_knee.mix(mapping_alevin_workflow.out.secondderiv_knee)
            ch_secondderiv_stats  = ch_secondderiv_stats.mix(mapping_alevin_workflow.out.secondderiv_stats)
            ch_secondderiv_cutoff = ch_secondderiv_cutoff.mix(mapping_alevin_workflow.out.secondderiv_cutoff)
        }

    emit:
        mapped_samplesheet           = ch_mapped_ss
        ref_gtf                      = ref_gtf_ch
        mapping_files                = ch_mapping_files
        starsolo_bam                 = ch_starsolo_bam
        star_solodir                 = ch_star_solodir
        starsolo_genefull50_raw      = ch_starsolo_genefull50_raw
        starsolo_genefull50_filtered = ch_starsolo_genefull50_filtered
        starsolo_velocyto_raw        = ch_starsolo_velocyto_raw
        starsolo_velocyto_filtered   = ch_starsolo_velocyto_filtered
        secondderiv_knee             = ch_secondderiv_knee
        secondderiv_stats            = ch_secondderiv_stats
        secondderiv_cutoff           = ch_secondderiv_cutoff
        saturation_imgs              = ch_sat_imgs
        saturation_residual_imgs     = ch_sat_res_imgs
        saturation_logs              = ch_sat_logs
        star_umipercell              = ch_star_umi
        star_log                     = ch_star_log
        star_final_log               = ch_star_final_log
        star_summaries               = ch_star_summaries
        star_cellreads               = ch_star_cellreads
        af_meta_info                 = ch_alevin_meta_info
        af_quant_json                = ch_alevin_quant_json
        af_cell_meta                 = ch_alevin_cell_meta
        af_mtx                       = ch_alevin_mtx
        af_filtered_mtx              = ch_alevin_filtered_mtx
        af_umipercell                = ch_alevin_umipercell
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
