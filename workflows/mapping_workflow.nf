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
include { mapping_starsolo_workflow as mapping_starsolo_geneext_input_workflow } from '../subworkflows/local/mapping/mapping_starsolo'
include { mapping_alevin_workflow                                           } from '../subworkflows/local/mapping/mapping_alevin'
include { mapping_alevin_workflow as mapping_alevin_geneext_workflow        } from '../subworkflows/local/mapping/mapping_alevin'
include { bam_inspection_workflow                                           } from '../subworkflows/local/post-processing/bam_inspection'
include { bam_inspection_workflow as bam_inspection_geneext_workflow        } from '../subworkflows/local/post-processing/bam_inspection'
include { bam_inspection_workflow as bam_inspection_geneext_input_workflow  } from '../subworkflows/local/post-processing/bam_inspection'
include { geneext_workflow                                                  } from '../subworkflows/local/mapping/geneext'

include { MERGE_REF_FASTA                                                   } from '../modules/local/custom/manipulate/merge_ref_fasta/main'
include { MERGE_REF_GTF                                                     } from '../modules/local/custom/manipulate/merge_ref_gtf/main'
include { MERGE_REF_GTF as MERGE_REF_GTF_GENEEXT                            } from '../modules/local/custom/manipulate/merge_ref_gtf/main'
include { FASTQC                                                            } from '../modules/local/tools/fastqc/main'
include { SUBSAMPLE_FASTQS                                                  } from '../modules/local/tools/seqtk/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN MAPPING
        Maps every analytical run lib/AnalyticalRuns.groovy plans for this combination of
        mapping_software, perform_geneext, geneext_downstream and run_method:

        - the standard runs, against the reference annotation;
        - the full-data STARsolo pass GeneExt reads, which is either the standard
          '_starsolo' run or, where that is not mapped, a pass of its own that only
          produces the BAM and is neither reported nor carried downstream;
        - the re-mappings against the GeneExt-extended annotation.
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
        def ch_mapped_ss                    = channel.empty()
        def ch_mapping_files                = channel.empty()
        def ch_starsolo_bam                 = channel.empty()
        def ch_star_solodir                 = channel.empty()
        def ch_starsolo_genefull50_raw      = channel.empty()
        def ch_starsolo_genefull50_filtered = channel.empty()
        def ch_starsolo_velocyto_raw        = channel.empty()
        def ch_starsolo_velocyto_filtered   = channel.empty()
        def ch_secondderiv_knee             = channel.empty()
        def ch_secondderiv_stats            = channel.empty()
        def ch_secondderiv_cutoff           = channel.empty()
        def ch_sat_imgs                     = channel.empty()
        def ch_sat_res_imgs                 = channel.empty()
        def ch_sat_logs                     = channel.empty()
        def ch_star_umi                     = channel.empty()
        def ch_star_log                     = channel.empty()
        def ch_star_final_log               = channel.empty()
        def ch_star_summaries               = channel.empty()
        def ch_star_cellreads               = channel.empty()
        def ch_alevin_meta_info             = channel.empty()
        def ch_alevin_quant_json            = channel.empty()
        def ch_alevin_cell_meta             = channel.empty()
        def ch_alevin_mtx                   = channel.empty()
        def ch_alevin_filtered_mtx          = channel.empty()
        def ch_alevin_umipercell            = channel.empty()
        def ch_featurecounts                = channel.empty()
        def ch_antisense                    = channel.empty()
        def ch_barcode_reads                = channel.empty()
        def ch_pavian_sankey                = channel.empty()
        def ch_kraken_report                = channel.empty()
        def ch_geneext_report               = channel.empty()
        def ch_geneext_log                  = channel.empty()
        def ch_versions                     = channel.empty()

        // The analytical runs to map, see lib/AnalyticalRuns.groovy
        def runs = AnalyticalRuns.plan(params)

        // Conditionally bypass MERGE_REF_GTF/FASTA when no additional features are provided
        def ref_gtf_ch
        if (params.ref_gtf_addfeature) {
            MERGE_REF_GTF(params.ref_gtf, channel.fromPath(params.ref_gtf_addfeature), '')
            ref_gtf_ch = MERGE_REF_GTF.out.gtf
            ch_versions = ch_versions.mix(MERGE_REF_GTF.out.versions)
        } else {
            ref_gtf_ch = channel.value(file(params.ref_gtf))
        }

        def ref_fasta_ch
        if (params.ref_fasta_addfeature) {
            MERGE_REF_FASTA(params.ref_fasta, channel.fromPath(params.ref_fasta_addfeature))
            ref_fasta_ch = MERGE_REF_FASTA.out.fasta
            ch_versions = ch_versions.mix(MERGE_REF_FASTA.out.versions)
        } else {
            ref_fasta_ch = channel.value(file(params.ref_fasta))
        }

        // The GeneExt reference, set once GeneExt has run, so both mappers re-map against
        // the same extended annotation
        def ref_gtf_geneext_ch   = null
        def ref_fasta_geneext_ch = null

        // Safe bc_whitelist: emit empty string when no whitelist is produced by preprocessing
        def bc_whitelist_safe = bc_whitelist.ifEmpty("")

        // Alevin takes the whitelist as a path input, which cannot be an empty string, so protocols running without a whitelist (e.g. MARS-seq) stage nothing instead
        def bc_whitelist_alevin = bc_whitelist_safe
            .map { wl -> wl?.toString()?.trim() ? wl : [] }

        // Quality Control
        FASTQC(ch_samples)
        ch_versions = ch_versions.mix(FASTQC.out.versions)

        // 'alevin_subsampled_starsolo' maps STARsolo on a subsample. It is drawn once and shared
        // by the standard and the GeneExt subsampled runs, so both map the same reads.
        def ch_subsampled = channel.empty()
        if (runs.star_standard == '_subsampled_starsolo' || runs.star_geneext == '_geneext_subsampled_starsolo') {
            SUBSAMPLE_FASTQS(apply_suffix(ch_samples, "_subsampled_starsolo"))
            ch_subsampled = SUBSAMPLE_FASTQS.out.subsampled_files
            ch_versions = ch_versions.mix(SUBSAMPLE_FASTQS.out.versions)
        }

        // The full-data STARsolo pass GeneExt reads, where it is not the standard run itself.
        // Only its BAM is used, so geneext_bam_only strips it down to the alignment.
        def star_index_standard = null
        if (runs.geneext_input == 'own_run') {
            mapping_starsolo_geneext_input_workflow(apply_suffix(ch_samples, AnalyticalRuns.GENEEXT_INPUT_SUFFIX),
                                                    bc_whitelist_safe, ref_gtf_ch, ref_fasta_ch, 'false',
                                                    params.geneext_bam_only as boolean, null)

            // Same reference, so a standard subsampled run does not build a second index
            star_index_standard = mapping_starsolo_geneext_input_workflow.out.star_index
            ch_versions = ch_versions.mix(mapping_starsolo_geneext_input_workflow.out.versions)

            // geneext_bam_only = false asks for this pass's own QC as well. It is published,
            // but not reported: the pass is not one of the analytical runs.
            if (!params.geneext_bam_only && params.star_generateBAM) {
                bam_inspection_geneext_input_workflow(mapping_starsolo_geneext_input_workflow.out.starsolo_bam, ref_gtf_ch,
                                                        mapping_starsolo_geneext_input_workflow.out.star_summaries,
                                                        mapping_starsolo_geneext_input_workflow.out.star_final_log,
                                                        mapping_starsolo_geneext_input_workflow.out.secondderiv_stats,
                                                        mapping_starsolo_geneext_input_workflow.out.starsolo_genefull50_filtered,
                                                        mapping_starsolo_geneext_input_workflow.out.star_cellreads)
                ch_versions = ch_versions.mix(bam_inspection_geneext_input_workflow.out.versions)
            }
        }

        // Standard STARsolo run, on the full data or on the subsample
        if (runs.star_standard) {
            def ch_star_input = runs.star_standard == '_subsampled_starsolo'
                ? ch_subsampled
                : apply_suffix(ch_samples, runs.star_standard)
            ch_mapped_ss = ch_mapped_ss.mix(apply_suffix(ch_samples, runs.star_standard))

            mapping_starsolo_workflow(ch_star_input, bc_whitelist_safe, ref_gtf_ch, ref_fasta_ch, 'false', false, star_index_standard)

            // Assign starsolo standard mapping outputs
            ch_mapping_files                = ch_mapping_files.mix(mapping_starsolo_workflow.out.mapping_files)
            ch_starsolo_bam                 = ch_starsolo_bam.mix(mapping_starsolo_workflow.out.starsolo_bam)
            ch_star_solodir                 = ch_star_solodir.mix(mapping_starsolo_workflow.out.star_solodir)
            ch_starsolo_genefull50_raw      = ch_starsolo_genefull50_raw.mix(mapping_starsolo_workflow.out.starsolo_genefull50_raw)
            ch_starsolo_genefull50_filtered = ch_starsolo_genefull50_filtered.mix(mapping_starsolo_workflow.out.starsolo_genefull50_filtered)
            ch_starsolo_velocyto_raw        = ch_starsolo_velocyto_raw.mix(mapping_starsolo_workflow.out.starsolo_velocyto_raw)
            ch_starsolo_velocyto_filtered   = ch_starsolo_velocyto_filtered.mix(mapping_starsolo_workflow.out.starsolo_velocyto_filtered)
            ch_secondderiv_knee             = ch_secondderiv_knee.mix(mapping_starsolo_workflow.out.secondderiv_knee)
            ch_secondderiv_stats            = ch_secondderiv_stats.mix(mapping_starsolo_workflow.out.secondderiv_stats)
            ch_secondderiv_cutoff           = ch_secondderiv_cutoff.mix(mapping_starsolo_workflow.out.secondderiv_cutoff)
            ch_star_umi                     = ch_star_umi.mix(mapping_starsolo_workflow.out.star_umipercell)
            ch_star_log                     = ch_star_log.mix(mapping_starsolo_workflow.out.star_log)
            ch_star_final_log               = ch_star_final_log.mix(mapping_starsolo_workflow.out.star_final_log)
            ch_star_summaries               = ch_star_summaries.mix(mapping_starsolo_workflow.out.star_summaries)
            ch_star_cellreads               = ch_star_cellreads.mix(mapping_starsolo_workflow.out.star_cellreads)
            ch_versions                     = ch_versions.mix(mapping_starsolo_workflow.out.versions)

            // Run BAM inspection workflow on STARsolo output
            if (params.star_generateBAM) {
                bam_inspection_workflow(mapping_starsolo_workflow.out.starsolo_bam, ref_gtf_ch,
                                            mapping_starsolo_workflow.out.star_summaries,
                                            mapping_starsolo_workflow.out.star_final_log,
                                            mapping_starsolo_workflow.out.secondderiv_stats,
                                            mapping_starsolo_workflow.out.starsolo_genefull50_filtered,
                                            mapping_starsolo_workflow.out.star_cellreads)

                ch_sat_imgs              =  bam_inspection_workflow.out.saturation_imgs
                ch_sat_res_imgs          =  bam_inspection_workflow.out.saturation_residual_imgs
                ch_sat_logs              =  bam_inspection_workflow.out.saturation_logs
                ch_featurecounts         =  bam_inspection_workflow.out.featurecount_txt
                ch_antisense             =  bam_inspection_workflow.out.antisense_txt
                ch_barcode_reads         =  bam_inspection_workflow.out.barcode_reads
                ch_pavian_sankey         =  bam_inspection_workflow.out.pavian_sankey
                ch_kraken_report         =  bam_inspection_workflow.out.kraken_report
                ch_versions              =  ch_versions.mix(bam_inspection_workflow.out.versions)
            }
        }

        // Extend the annotation with GeneExt, from the full-data STARsolo alignments
        if (runs.geneext_input) {

            def ch_geneext_input_bam = runs.geneext_input == 'reuse_standard'
                ? mapping_starsolo_workflow.out.starsolo_bam
                : mapping_starsolo_geneext_input_workflow.out.starsolo_bam

            // Collect ALL BAMs from all samples before running geneext
            geneext_workflow(ch_geneext_input_bam)

            // GeneExt's own run statistics, summarised in the dashboard
            ch_geneext_report = geneext_workflow.out.report
            ch_geneext_log    = geneext_workflow.out.geneext_log
            ch_versions       = ch_versions.mix(geneext_workflow.out.versions)

            if (runs.star_geneext || runs.alevin_geneext) {

                // Same conditional bypass for geneext GTF
                if (params.ref_gtf_addfeature) {
                    MERGE_REF_GTF_GENEEXT(geneext_workflow.out.ref_gtf, channel.fromPath(params.ref_gtf_addfeature), '_geneext')
                    ref_gtf_geneext_ch = MERGE_REF_GTF_GENEEXT.out.gtf
                } else {
                    // Geneext always extends from the geneext output, no bypass possible here
                    MERGE_REF_GTF_GENEEXT(geneext_workflow.out.ref_gtf, channel.value([]), '_geneext')
                    ref_gtf_geneext_ch = MERGE_REF_GTF_GENEEXT.out.gtf
                }
                ch_versions = ch_versions.mix(MERGE_REF_GTF_GENEEXT.out.versions)

                // GeneExt only extends the annotation, so the genome is the standard one
                ref_fasta_geneext_ch = ref_fasta_ch
            }
        }

        // STARsolo re-mapped against the GeneExt annotation, on the same reads as the standard
        // STARsolo run of this mapping_software would use
        if (runs.star_geneext) {
            def ch_geneext = runs.star_geneext == '_geneext_subsampled_starsolo'
                ? ch_subsampled.map { row ->
                    def new_meta = row[0].clone()
                    new_meta.id = row[0].base_id + runs.star_geneext
                    [new_meta] + row.drop(1)
                }
                : apply_suffix(ch_samples, runs.star_geneext)
            ch_mapped_ss = ch_mapped_ss.mix(apply_suffix(ch_samples, runs.star_geneext))

            mapping_starsolo_geneext_workflow(ch_geneext, bc_whitelist_safe, ref_gtf_geneext_ch, ref_fasta_geneext_ch, 'true', false, null)

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
            ch_versions                     = ch_versions.mix(mapping_starsolo_geneext_workflow.out.versions)

            // Run BAM inspection workflow on geneext STARsolo output
            if (params.star_generateBAM) {
                bam_inspection_geneext_workflow(mapping_starsolo_geneext_workflow.out.starsolo_bam,
                                            ref_gtf_geneext_ch,
                                            mapping_starsolo_geneext_workflow.out.star_summaries,
                                            mapping_starsolo_geneext_workflow.out.star_final_log,
                                            mapping_starsolo_geneext_workflow.out.secondderiv_stats,
                                            mapping_starsolo_geneext_workflow.out.starsolo_genefull50_filtered,
                                            mapping_starsolo_geneext_workflow.out.star_cellreads)

                ch_sat_imgs                     = ch_sat_imgs.mix(bam_inspection_geneext_workflow.out.saturation_imgs)
                ch_sat_res_imgs                 = ch_sat_res_imgs.mix(bam_inspection_geneext_workflow.out.saturation_residual_imgs)
                ch_sat_logs                     = ch_sat_logs.mix(bam_inspection_geneext_workflow.out.saturation_logs)
                ch_featurecounts                = ch_featurecounts.mix(bam_inspection_geneext_workflow.out.featurecount_txt)
                ch_antisense                    = ch_antisense.mix(bam_inspection_geneext_workflow.out.antisense_txt)
                ch_barcode_reads                = ch_barcode_reads.mix(bam_inspection_geneext_workflow.out.barcode_reads)
                ch_pavian_sankey                = ch_pavian_sankey.mix(bam_inspection_geneext_workflow.out.pavian_sankey)
                ch_kraken_report                = ch_kraken_report.mix(bam_inspection_geneext_workflow.out.kraken_report)
                ch_versions                     = ch_versions.mix(bam_inspection_geneext_workflow.out.versions)
            }
        }

        // Standard alevin-fry run
        if (runs.alevin_standard) {
            def ch_alevin = apply_suffix(ch_samples, "_alevinfry")
            ch_mapped_ss = ch_mapped_ss.mix(ch_alevin)
            // Same reference as the standard STARsolo run, so the two mappers are compared like for like
            mapping_alevin_workflow(ch_alevin, bc_whitelist_alevin, ref_gtf_ch, ref_fasta_ch, 'false')

            ch_mapping_files       = ch_mapping_files.mix(mapping_alevin_workflow.out.mapping_files)
            ch_alevin_meta_info    = ch_alevin_meta_info.mix(mapping_alevin_workflow.out.af_meta_info)
            ch_alevin_quant_json   = ch_alevin_quant_json.mix(mapping_alevin_workflow.out.af_quant_json)
            ch_alevin_cell_meta    = ch_alevin_cell_meta.mix(mapping_alevin_workflow.out.af_cell_meta)
            ch_alevin_mtx          = ch_alevin_mtx.mix(mapping_alevin_workflow.out.af_mtx)
            ch_alevin_filtered_mtx = ch_alevin_filtered_mtx.mix(mapping_alevin_workflow.out.af_filtered_mtx)
            ch_alevin_umipercell   = ch_alevin_umipercell.mix(mapping_alevin_workflow.out.af_umipercell)

            // Both mappers emit the same second-derivative artefacts, so the reporting channels carry them together
            ch_secondderiv_knee   = ch_secondderiv_knee.mix(mapping_alevin_workflow.out.secondderiv_knee)
            ch_secondderiv_stats  = ch_secondderiv_stats.mix(mapping_alevin_workflow.out.secondderiv_stats)
            ch_secondderiv_cutoff = ch_secondderiv_cutoff.mix(mapping_alevin_workflow.out.secondderiv_cutoff)

            ch_versions = ch_versions.mix(mapping_alevin_workflow.out.versions)
        }

        // alevin-fry re-mapped against the GeneExt annotation, the counterpart of the GeneExt STARsolo run
        if (runs.alevin_geneext) {
            def ch_alevin_geneext = apply_suffix(ch_samples, "_geneext_alevinfry")
            ch_mapped_ss = ch_mapped_ss.mix(ch_alevin_geneext)
            mapping_alevin_geneext_workflow(ch_alevin_geneext, bc_whitelist_alevin, ref_gtf_geneext_ch, ref_fasta_geneext_ch, 'true')

            ch_mapping_files       = ch_mapping_files.mix(mapping_alevin_geneext_workflow.out.mapping_files)
            ch_alevin_meta_info    = ch_alevin_meta_info.mix(mapping_alevin_geneext_workflow.out.af_meta_info)
            ch_alevin_quant_json   = ch_alevin_quant_json.mix(mapping_alevin_geneext_workflow.out.af_quant_json)
            ch_alevin_cell_meta    = ch_alevin_cell_meta.mix(mapping_alevin_geneext_workflow.out.af_cell_meta)
            ch_alevin_mtx          = ch_alevin_mtx.mix(mapping_alevin_geneext_workflow.out.af_mtx)
            ch_alevin_filtered_mtx = ch_alevin_filtered_mtx.mix(mapping_alevin_geneext_workflow.out.af_filtered_mtx)
            ch_alevin_umipercell   = ch_alevin_umipercell.mix(mapping_alevin_geneext_workflow.out.af_umipercell)
            ch_secondderiv_knee    = ch_secondderiv_knee.mix(mapping_alevin_geneext_workflow.out.secondderiv_knee)
            ch_secondderiv_stats   = ch_secondderiv_stats.mix(mapping_alevin_geneext_workflow.out.secondderiv_stats)
            ch_secondderiv_cutoff  = ch_secondderiv_cutoff.mix(mapping_alevin_geneext_workflow.out.secondderiv_cutoff)
            ch_versions            = ch_versions.mix(mapping_alevin_geneext_workflow.out.versions)
        }

    emit:
        fastqc_results               = FASTQC.out.fastqc_results
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
        antisense_txt                = ch_antisense
        barcode_reads                = ch_barcode_reads
        pavian_sankey                = ch_pavian_sankey
        kraken_report                = ch_kraken_report
        geneext_report               = ch_geneext_report
        geneext_log                  = ch_geneext_log
        versions                     = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
