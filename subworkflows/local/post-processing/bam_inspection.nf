//
// Subworkflow with functionality specific to the workflow 'mapping_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_INDEX                                    } from '../../../modules/local/tools/samtools/samtools_index/main'
include { SATURATION_TABLE                                  } from '../../../modules/local/tools/10x_saturate/saturation_table/main'
include { SATURATION_PLOT                                   } from '../../../modules/local/tools/10x_saturate/plot_curve/main'
include { SAMTOOLS_VIEW_MAPPED                              } from '../../../modules/local/tools/samtools/samtools_view_mapped/main'
include { SAMTOOLS_VIEW_UNMAPPED                            } from '../../../modules/local/tools/samtools/samtools_view_unmapped/main'
include { CALC_MT_RRNA                                      } from '../../../modules/local/tools/featurecounts/main'
include { KRAKEN_CREATE_DB                                  } from '../../../modules/local/tools/kraken/kraken_create_db/main'
include { KRAKEN                                            } from '../../../modules/local/tools/kraken/kraken_classify/main'
include { PAVIAN                                            } from '../../../modules/local/tools/pavian/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN ALEVIN-FRY MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow bam_inspection_workflow {
    take:
        bam_file
        ref_gtf
        summary_csv
        log_final_file

    main:
        // Initialize reporting channels
        def ch_sat_imgs                 = Channel.empty()
        def ch_sat_res_imgs             = Channel.empty()
        def ch_sat_logs                 = Channel.empty()
        def ch_featurecounts            = Channel.empty()
        def ch_pavian_sankey            = Channel.empty()

        SAMTOOLS_INDEX(bam_file)

        // Calculate saturation curve if perform_10x_saturate is true
        if (params.perform_10x_saturate) {
            SAMTOOLS_VIEW_MAPPED(bam_file)

            // Join channels on sample ID before 10x_saturate
            SAMTOOLS_VIEW_MAPPED.out.filtered_mapped_bam
                .join(summary_csv)
                .join(log_final_file)
                .join(SAMTOOLS_VIEW_MAPPED.out.filtered_mapped_bai)
                .join(SAMTOOLS_VIEW_MAPPED.out.mapreads)
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
            bam_file
                .join(SAMTOOLS_INDEX.out.bam_index)
                .multiMap { meta, bam_file, bai ->
                    bam_ch: [meta, bam_file]
                    bai_ch:  [meta, bai]
                }
                .set { ch_fc_inputs }

            // Run featureCounts to calculate mtDNA and rRNA percentages and capture output
            CALC_MT_RRNA(ch_fc_inputs.bam_ch, ch_fc_inputs.bai_ch, ref_gtf.first())
            ch_featurecounts = CALC_MT_RRNA.out
        }

        // Inspecting unmapped reads using Kraken2
        if (params.perform_kraken) {
            // Extract unmapped reads
            SAMTOOLS_VIEW_UNMAPPED(bam_file)

            // Perform kraken tools plus visualization with Pavian
            KRAKEN_CREATE_DB()
            KRAKEN(KRAKEN_CREATE_DB.out.db_path_file, SAMTOOLS_VIEW_UNMAPPED.out.filtered_unmapped_fasta)
            PAVIAN(KRAKEN.out)
            ch_pavian_sankey = PAVIAN.out
        }

    emit:
        saturation_imgs                 = ch_sat_imgs
        saturation_residual_imgs        = ch_sat_res_imgs
        saturation_logs                 = ch_sat_logs
        featurecount_txt                = ch_featurecounts
        pavian_sankey                   = ch_pavian_sankey
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
