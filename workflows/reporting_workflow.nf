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
        run_config
        star_logs
        star_summaries
        star_full_logs
        starsolo_bam
        star_solodir
        saturation_logs
        cell_stats
        af_meta_info
        af_quant_json
        af_cell_meta
        sankey_files
        saturation_imgs
        residuals_imgs
        knee_files
        mt_rrna_metrics

    main:
        // Join BAM, SoloDir, and Logs before per-cell metrics
        starsolo_bam
            .join(star_solodir)
            .join(star_logs)
            .multiMap { meta, bam, solodir, logs ->
                bam_ch:     [meta, bam]
                solodir_ch: [meta, solodir]
                logs_ch:    [meta, logs]
            }
            .set { ch_percell_inputs }

        // Only create per-cell plots if params.mt_contig is set by user
        if (params.mt_contig != "chrM M MT") {
            PERCELL_METRICS(
                ch_percell_inputs.bam_ch,
                ch_percell_inputs.solodir_ch,
                ch_percell_inputs.logs_ch
            )
            percell_json = PERCELL_METRICS.out.percell_json
        } else {
            percell_json = Channel.empty()
        }

        // Join channels that need renaming by meta ID
        ch_to_rename = samplesheet
            .map { row -> [ row[0] ] }
            .join(star_summaries, remainder: true)
            .join(cell_stats, remainder: true)
            .join(knee_files, remainder: true)
            .join(af_meta_info, remainder: true)
            .join(af_quant_json, remainder: true)
            .join(af_cell_meta, remainder: true)
            .map { row ->
                def meta = row[0]
                // Extract all elements after 'meta', flatten them, and filter out nulls
                def input_files = row[1..-1].flatten().findAll { it != null }
                [ meta, input_files ]
            }

        // Rename STARsolo files
        PREPARE_DASHBOARD_INPUTS(ch_to_rename)

        // Run Dashboard Generation
        GENERATE_DASHBOARD(
            samplesheet_file,
            run_config,
            star_logs.map{ it[1] }.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.summary.collect().ifEmpty([]),
            star_full_logs.map{ it[1] }.collect().ifEmpty([]),
            saturation_logs.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.cell_stats.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.af_meta_info.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.af_quant_json.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.af_cell_meta.collect().ifEmpty([]),
            sankey_files.collect().ifEmpty([]),
            saturation_imgs.collect().ifEmpty([]),
            residuals_imgs.collect().ifEmpty([]),
            PREPARE_DASHBOARD_INPUTS.out.knee_files.collect().ifEmpty([]),
            mt_rrna_metrics.collect().ifEmpty([]),
            percell_json.collect().ifEmpty([])
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
        multiqc_report  = MULTIQC.out
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
