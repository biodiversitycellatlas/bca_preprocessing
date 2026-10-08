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
include { CALC_READ_METRICS                                 } from '../../../modules/local/custom/read_metrics/main'
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
        secondderiv_stats
        filtered_matrix     // tuple(meta, filtered matrix dir); restricts the called-cell metrics
        cellreads_stats     // tuple(meta, GeneFull_Ex50pAS/CellReads.stats); source of the antisense metrics

    main:
        // Initialize reporting channels
        def ch_sat_imgs                 = Channel.empty()
        def ch_sat_res_imgs             = Channel.empty()
        def ch_sat_logs                 = Channel.empty()
        def ch_pavian_sankey            = Channel.empty()
        def ch_kraken_report            = Channel.empty()

        SAMTOOLS_INDEX(bam_file)
        def ch_versions = SAMTOOLS_INDEX.out.versions

        // Calculate saturation curve if perform_10x_saturate is true
        if (params.perform_10x_saturate) {
            SAMTOOLS_VIEW_MAPPED(bam_file)

            // Join channels on sample ID before 10x_saturate
            SAMTOOLS_VIEW_MAPPED.out.filtered_mapped_bam
                .join(summary_csv)
                .join(log_final_file)
                .join(SAMTOOLS_VIEW_MAPPED.out.filtered_mapped_bai)
                .join(SAMTOOLS_VIEW_MAPPED.out.mapreads)
                .join(secondderiv_stats, remainder: true)
                .filter { row -> row.size() == 7 && row[1] != null }
                .multiMap { meta, bam, summary, log_final, bai, mapreads, sd_stats ->
                    bam_ch:         [meta, bam]
                    summary_ch:     [meta, summary]
                    log_final_ch:   [meta, log_final]
                    bai_ch:         [meta, bai]
                    reads_ch:       [meta, mapreads]
                    sd_stats_ch:    [meta, sd_stats ?: []]
                }
                .set { ch_saturation_inputs }

            SATURATION_TABLE(
                ch_saturation_inputs.bam_ch,
                ch_saturation_inputs.summary_ch,
                ch_saturation_inputs.log_final_ch,
                ch_saturation_inputs.bai_ch,
                ch_saturation_inputs.reads_ch,
                ch_saturation_inputs.sd_stats_ch
            )

            SATURATION_PLOT(SATURATION_TABLE.out.saturation_table)

            // Capture saturation plot outputs
            ch_sat_imgs     = SATURATION_PLOT.out.img_saturation
            ch_sat_res_imgs = SATURATION_PLOT.out.img_residuals
            ch_sat_logs     = SATURATION_PLOT.out.logs

            ch_versions = ch_versions.mix(
                SAMTOOLS_VIEW_MAPPED.out.versions,
                SATURATION_TABLE.out.versions,
                SATURATION_PLOT.out.versions
            )
        }

        // Percentages of mtDNA and rRNA reads, per library, per called cell and per barcode,
        // in one pass over the BAM, and the antisense share of gene reads from STARsolo's
        // CellReads.stats (stranded runs only). Always run: PERCELL_METRICS reads its
        // per-barcode counts from here instead of scanning the BAM a second time.
        //
        // Join each BAM with, where they exist, the filtered matrix (called cells) and
        // CellReads.stats. remainder keeps BAMs without them (bam_only runs), but can
        // also emit rows for a matrix without a BAM; those are dropped here
        bam_file
            .join(filtered_matrix, remainder: true)
            .join(cellreads_stats, remainder: true)
            .filter { row -> row.size() == 4 && row[1] != null }
            .map { meta, bam, filtered, cellreads -> [meta, bam, filtered ?: [], cellreads ?: []] }
            .set { ch_metrics_inputs }

        // The added features are the rRNA reference, so every read on their contigs
        // counts as rRNA. Read from params rather than taken as a subworkflow input,
        // since ref_gtf already arrives merged with this file and the module needs
        // the two apart to tell the added contigs from the rest.
        def ch_rrna_gtf = params.ref_gtf_addfeature
            ? Channel.value(file(params.ref_gtf_addfeature))
            : Channel.value([])

        CALC_READ_METRICS(ch_metrics_inputs, ref_gtf.first(), ch_rrna_gtf)
        ch_versions = ch_versions.mix(CALC_READ_METRICS.out.versions)

        // Inspecting unmapped reads using Kraken2
        if (params.perform_kraken) {

            // Extract unmapped reads
            SAMTOOLS_VIEW_UNMAPPED(bam_file)

            // Perform kraken tools plus visualization with Pavian
            KRAKEN_CREATE_DB()
            KRAKEN(KRAKEN_CREATE_DB.out.db_path_file, SAMTOOLS_VIEW_UNMAPPED.out.filtered_unmapped_fasta)
            PAVIAN(KRAKEN.out.k2report)
            ch_pavian_sankey = PAVIAN.out.sankey
            ch_kraken_report = KRAKEN.out.k2report

            ch_versions = ch_versions.mix(
                SAMTOOLS_VIEW_UNMAPPED.out.versions,
                KRAKEN_CREATE_DB.out.versions,
                KRAKEN.out.versions,
                PAVIAN.out.versions
            )
        }

    emit:
        saturation_imgs                 = ch_sat_imgs
        saturation_residual_imgs        = ch_sat_res_imgs
        saturation_logs                 = ch_sat_logs
        featurecount_txt                = CALC_READ_METRICS.out.mt_rrna_metrics
        antisense_txt                   = CALC_READ_METRICS.out.antisense_metrics
        barcode_reads                   = CALC_READ_METRICS.out.barcode_reads
        pavian_sankey                   = ch_pavian_sankey
        kraken_report                   = ch_kraken_report
        versions                        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
