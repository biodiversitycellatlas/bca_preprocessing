//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MAPPING_STATS             } from '../modules/local/custom/dashboard/mapping_stats/main'
include { MULTIQC                   } from '../modules/local/tools/multiqc/main'
include { PREPARE_DASHBOARD_INPUTS  } from '../modules/local/custom/dashboard/prepare_inputs/main'
include { GENERATE_DASHBOARD        } from '../modules/local/custom/dashboard/html_report/main'
include { PERCELL_METRICS           } from '../modules/local/custom/dashboard/per_cell_metrics/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN FILTERING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow reporting_workflow {
    take:
        samplesheet
        samplesheet_file
        ref_gtf
        run_config
        star_logs
        star_summaries
        star_full_logs
        starsolo_bam
        star_solodir
        star_filtered_mtx
        saturation_logs
        cell_stats
        af_meta_info
        af_quant_json
        af_cell_meta
        af_mtx
        af_umipercell
        sankey_files
        saturation_imgs
        residuals_imgs
        knee_files
        mt_rrna_metrics
        secondderiv_knee
        secondderiv_stats
        secondderiv_cutoff
        cs_ambient_hist_plot
        cs_umap_comparison_plot
        cs_top_genes
        geneext_report
        geneext_log

    main:
        // Join BAM, SoloDir, and Logs before per-cell metrics. The cutoff and the
        // filtered matrix are both optional: the module prefers the filtered matrix's
        // barcodes, falls back to the cutoff, and then to STARsolo's nUMImin.
        starsolo_bam
            .join(star_solodir)
            .join(star_logs)
            .join(secondderiv_cutoff, remainder: true)
            .join(star_filtered_mtx, remainder: true)
            // remainder keeps samples without those inputs, but can also emit rows for
            // samples that have them and no BAM; those are dropped here
            .filter { row -> row.size() == 6 && row[1] != null }
            .multiMap { meta, bam, solodir, logs, cutoff, filtered ->
                bam_ch:      [meta, bam]
                solodir_ch:  [meta, solodir]
                logs_ch:     [meta, logs]
                cutoff_ch:   [meta, cutoff ?: []]
                filtered_ch: [meta, filtered ?: []]
            }
            .set { ch_percell_inputs }

        // Run per-cell metrics on starsolo outputs.
        PERCELL_METRICS(
            ch_percell_inputs.bam_ch,
            ch_percell_inputs.solodir_ch,
            ch_percell_inputs.logs_ch,
            ch_percell_inputs.cutoff_ch,
            ch_percell_inputs.filtered_ch,
            ref_gtf.first()
        )
        percell_json = PERCELL_METRICS.out.percell_json

        // Build anchor from whichever mapping software was run
        def ch_anchor = star_summaries.mix(af_meta_info)

        // quant.json's num_genes counts USA matrix columns (3 per gene)
        def ch_af_mat_cols = af_mtx.map { meta, dir ->
            [ meta, file("${dir}/quants_mat_cols.txt") ]
        }

        // Join channels that need renaming by meta ID
        ch_to_rename = ch_anchor
            .join(cell_stats,     remainder: true)
            .join(knee_files,     remainder: true)
            .join(af_meta_info,   remainder: true)
            .join(af_quant_json,  remainder: true)
            .join(af_cell_meta,   remainder: true)
            .join(ch_af_mat_cols, remainder: true)
            .join(af_umipercell,  remainder: true)
            .map { row ->
                def meta        = row[0]
                def input_files = row.drop(1).findAll { item ->
                    item != null && !(item instanceof java.util.Map)
                }
                [ meta, input_files ]
            }

        // Rename STARsolo files
        PREPARE_DASHBOARD_INPUTS(ch_to_rename)

        // Extract the exact analytical mappings that were successfully processed
        def ch_analytical_manifest = star_logs
            .map { meta, log -> "${meta.id},${meta.base_id ?: meta.id},starsolo" }
            .mix(
                af_meta_info.map { meta, info -> "${meta.id},${meta.base_id ?: meta.id},alevin" }
            )
            .unique()
            .collectFile(
                name: 'analytical_samples.csv',
                newLine: true,
                seed: "analytical_id,base_id,source\n"
            )

        // Run Dashboard Generation
        GENERATE_DASHBOARD(
            samplesheet_file,
            run_config,
            ch_analytical_manifest,
            star_logs.map{ it[1] }.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.summary.collect().ifEmpty([]),
            star_full_logs.map{ it[1] }.collect().ifEmpty([]),
            saturation_logs.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.cell_stats.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.af_meta_info.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.af_quant_json.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.af_cell_meta.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.af_mat_cols.collect().ifEmpty([]),
            sankey_files.collect().ifEmpty([]),
            saturation_imgs.collect().ifEmpty([]),
            residuals_imgs.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.knee_files.collect().ifEmpty([]),
            mt_rrna_metrics.collect().ifEmpty([]),
            secondderiv_knee.map { it[1] }.collect().ifEmpty([]),
            secondderiv_stats.map { it[1] }.collect().ifEmpty([]),
            percell_json.collect().ifEmpty([]),
            cs_ambient_hist_plot.collect().ifEmpty([]),
            cs_umap_comparison_plot.collect().ifEmpty([]),
            cs_top_genes.collect().ifEmpty([]),
            geneext_report.collect().ifEmpty([]),
            geneext_log.collect().ifEmpty([])
        )

        // Trigger Mapping Stats & MultiQC
        mapping_stats_trigger = star_logs.collect()
            .mix(
                saturation_logs.collect(),
                mt_rrna_metrics.collect(),
                percell_json.collect()
            )
            .collect()
            .map { it -> true }

        MAPPING_STATS(mapping_stats_trigger)

        ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        MULTIQC(
            mapping_stats_trigger,
            ch_multiqc_config
        )

    emit:
        dashboard_html  = GENERATE_DASHBOARD.out.html
        multiqc_report  = MULTIQC.out.report
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
