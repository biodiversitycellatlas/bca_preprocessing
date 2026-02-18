//
// Subworkflow with functionality specific to the workflow 'mapping_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { STARSOLO_INDEX as STARSOLO_INDEX                  } from '../../../modules/local/tools/star/starsolo_genome_generate/main'
include { STARSOLO_INDEX as STARSOLO_INDEX_GENEEXT          } from '../../../modules/local/tools/star/starsolo_genome_generate/main'
include { STARSOLO_ALIGN as STARSOLO_ALIGN                  } from '../../../modules/local/tools/star/starsolo_align/main'
include { STARSOLO_ALIGN as STARSOLO_ALIGN_GENEEXT          } from '../../../modules/local/tools/star/starsolo_align/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX                  } from '../../../modules/local/tools/samtools/samtools_index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_GENEEXT          } from '../../../modules/local/tools/samtools/samtools_index/main'
include { SAMTOOLS_VIEW                                     } from '../../../modules/local/tools/samtools/samtools_view/main'
include { SATURATION_TABLE                                  } from '../../../modules/local/tools/10x_saturate/saturation_table/main'
include { SATURATION_PLOT                                   } from '../../../modules/local/tools/10x_saturate/plot_curve/main'
include { CALC_MT_RRNA as CALC_MT_RRNA                      } from '../../../modules/local/tools/featurecounts/main'
include { CALC_MT_RRNA as CALC_MT_RRNA_GENEEXT              } from '../../../modules/local/tools/featurecounts/main'
include { GENE_EXT                                          } from '../../../modules/local/tools/geneext/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN STARSOLO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow mapping_starsolo_workflow {
    take:
        data_output
        bc_whitelist

    main:
        // Initialize reporting channels
        def ch_mapping_files    = Channel.empty()
        def ch_starsolo_bam     = Channel.empty()
        def ch_star_solodir     = Channel.empty()
        def ch_genefull50_raw_dir = Channel.empty()
        def ch_sat_imgs         = Channel.empty()
        def ch_sat_res_imgs     = Channel.empty()
        def ch_sat_logs         = Channel.empty()
        def ch_star_umi         = Channel.empty()
        def ch_star_log         = Channel.empty()
        def ch_star_final_log   = Channel.empty()
        def ch_star_summaries   = Channel.empty()
        def ch_star_cellreads   = Channel.empty()
        def ch_featurecounts    = Channel.empty()

        // Get the first cDNA fastq file from the samplesheet to use for STAR index generation
        def ch_first_cDNA = data_output
            .map { meta, fastq_cDNA, fastq_BC_UMI, empty_list, samplesheet -> fastq_cDNA }
            .first()

        // Check if star index is provided, if not create it
        def star_index_ch
        if (params.star_index && file(params.star_index).exists()) {
            star_index_ch = Channel.value(file(params.star_index))
        } else {
            STARSOLO_INDEX(params.ref_gtf, ch_first_cDNA)
            star_index_ch = STARSOLO_INDEX.out
        }

        // Run STARsolo alignment
        STARSOLO_ALIGN(data_output, bc_whitelist, star_index_ch)
        ch_mapping_files = STARSOLO_ALIGN.out.starsolo_files

        // Capture STARsolo outputs
        ch_star_umi         = STARSOLO_ALIGN.out.umi_per_cell
        ch_star_log         = STARSOLO_ALIGN.out.log_file
        ch_star_final_log   = STARSOLO_ALIGN.out.log_final_file
        ch_star_summaries   = STARSOLO_ALIGN.out.summary_csv
        ch_star_cellreads   = STARSOLO_ALIGN.out.cellreads_stats
        ch_star_solodir     = STARSOLO_ALIGN.out.star_solodir
        ch_genefull50_raw_dir = STARSOLO_ALIGN.out.genefull50_raw_dir


        // Optional saturation analysis, feature inspection and gene extension if BAM is generated
        if (params.star_generateBAM) {

            // Capture the BAM output if generated and index
            ch_starsolo_bam = STARSOLO_ALIGN.out.bam_file
            SAMTOOLS_INDEX(ch_starsolo_bam)

            // Calculate saturation curve if perform_10x_saturate is true
            if (params.perform_10x_saturate) {
                SAMTOOLS_VIEW(STARSOLO_ALIGN.out.starsolo_files)

                // Join channels on sample ID before 10x_saturate
                SAMTOOLS_VIEW.out.filtered_bam
                    .join(STARSOLO_ALIGN.out.summary_csv)
                    .join(STARSOLO_ALIGN.out.log_final_file)
                    .join(SAMTOOLS_VIEW.out.filtered_bam_index)
                    .join(SAMTOOLS_VIEW.out.mapreads)
                    .multiMap { meta, bam, summary, log_final, bai, mapreads ->
                        bam_ch:         [meta, bam]
                        summary_ch:     [meta, summary]
                        log_final_ch:   [meta, log_final]
                        bai_ch:         [meta, bai]
                        reads_ch:       [meta, mapreads]
                    }
                    .set { ch_saturation_inputs }

                SATURATION_TABLE(
                    ch_saturation_inputs.bam_ch,
                    ch_saturation_inputs.summary_ch,
                    ch_saturation_inputs.log_final_ch,
                    ch_saturation_inputs.bai_ch,
                    ch_saturation_inputs.reads_ch
                )

                SATURATION_PLOT(SATURATION_TABLE.out)

                // Capture saturation plot outputs
                ch_sat_imgs     = SATURATION_PLOT.out.img_saturation
                ch_sat_res_imgs = SATURATION_PLOT.out.img_residuals
                ch_sat_logs     = SATURATION_PLOT.out.logs
            }

            // Calculate percentages mitochondrial DNA and ribosomal RNA
            if (params.perform_featurecounts) {
                // Join STARsolo files with samtools index
                STARSOLO_ALIGN.out.starsolo_files
                    .join(SAMTOOLS_INDEX.out.bam_index)
                    .multiMap { meta, star_files, bai ->
                        star_ch: [meta, star_files]
                        bai_ch:  [meta, bai]
                    }
                    .set { ch_fc_inputs }

                // Run featureCounts to calculate mtDNA and rRNA percentages and capture output
                CALC_MT_RRNA(ch_fc_inputs.star_ch, ch_fc_inputs.bai_ch)
                ch_featurecounts = CALC_MT_RRNA.out
            }

            // Gene Extension
            if (params.perform_geneext || params.run_method == "geneext_only") {

                // Join inputs for GENE_EXT
                STARSOLO_ALIGN.out.starsolo_files
                    .join(SAMTOOLS_INDEX.out.bam_index)
                    .multiMap { meta, star_files, bai ->
                        star_ch: [meta, star_files]
                        bai_ch:  [meta, bai]
                    }
                    .set { ch_geneext_inputs }

                // Run gene extension using GeneExt
                GENE_EXT(ch_geneext_inputs.star_ch, ch_geneext_inputs.bai_ch)

                // Remap STARsolo with extended GTF if run_method is not "geneext_only"
                if (params.run_method != "geneext_only") {

                    // Create STAR index with extended GTF
                    STARSOLO_INDEX_GENEEXT(GENE_EXT.out, ch_first_cDNA)

                    // Join original data with the new sample-specific index
                    data_output
                        .join(STARSOLO_INDEX_GENEEXT.out)
                        .multiMap { meta, f_cdna, f_umi, empty, sheet, index ->
                            data_ch:  [meta, f_cdna, f_umi, empty, sheet]
                            index_ch: [meta, index]
                        }
                        .set { ch_remap_inputs }

                    // Remap with STARsolo using the extended GTF
                    STARSOLO_ALIGN_GENEEXT(ch_remap_inputs.data_ch, bc_whitelist, ch_remap_inputs.index_ch)
                    SAMTOOLS_INDEX_GENEEXT(STARSOLO_ALIGN_GENEEXT.out.bam_file)

                    // Capture remapped STARsolo outputs
                    ch_mapping_files = ch_mapping_files.mix(STARSOLO_ALIGN_GENEEXT.out.starsolo_files)
                    ch_star_umi       = ch_star_umi.mix(STARSOLO_ALIGN_GENEEXT.out.umi_per_cell)
                    ch_star_log       = ch_star_log.mix(STARSOLO_ALIGN_GENEEXT.out.log_file)
                    ch_star_final_log = ch_star_final_log.mix(STARSOLO_ALIGN_GENEEXT.out.log_final_file)
                    ch_star_summaries = ch_star_summaries.mix(STARSOLO_ALIGN_GENEEXT.out.summary_csv)
                    ch_star_cellreads = ch_star_cellreads.mix(STARSOLO_ALIGN_GENEEXT.out.cellreads_stats)
                }
            }
        }

    emit:
        mapping_files            = ch_mapping_files
        starsolo_bam             = ch_starsolo_bam
        star_solodir             = ch_star_solodir
        starsolo_genefull50_raw  = ch_genefull50_raw_dir
        saturation_imgs          = ch_sat_imgs
        saturation_residual_imgs = ch_sat_res_imgs
        saturation_logs          = ch_sat_logs
        star_umipercell          = ch_star_umi
        star_log                 = ch_star_log
        star_final_log           = ch_star_final_log
        star_summaries           = ch_star_summaries
        star_cellreads           = ch_star_cellreads
        featurecount_txt         = ch_featurecounts
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
