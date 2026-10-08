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
        barcode_reads
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
        antisense_metrics
        secondderiv_knee
        secondderiv_stats
        secondderiv_cutoff
        cs_ambient_hist_plot
        cs_umap_comparison_plot
        cs_top_genes
        geneext_report
        geneext_log
        fastqc_results
        kraken_report
        starsolo_raw_mtx
        vendor_results

    main:
        // Declare channels
        def ch_star_logs = star_logs

        // Join the per-barcode read counts (CALC_READ_METRICS, one per STARsolo BAM),
        // SoloDir, and Logs before per-cell metrics. The cutoff and the filtered matrix
        // are both optional: the module prefers the filtered matrix's barcodes, falls
        // back to the cutoff, and then to STARsolo's nUMImin.
        barcode_reads
            .join(star_solodir)
            .join(ch_star_logs)
            .join(secondderiv_cutoff, remainder: true)
            .join(star_filtered_mtx, remainder: true)
            // remainder keeps samples without those inputs, but can also emit rows for
            // samples that have them and no read counts; those are dropped here
            .filter { row -> row.size() == 6 && row[1] != null }
            .multiMap { meta, counts, solodir, logs, cutoff, filtered ->
                counts_ch:   [meta, counts]
                solodir_ch:  [meta, solodir]
                logs_ch:     [meta, logs]
                cutoff_ch:   [meta, cutoff ?: []]
                filtered_ch: [meta, filtered ?: []]
            }
            .set { ch_percell_inputs }

        // Run per-cell metrics on starsolo outputs.
        PERCELL_METRICS(
            ch_percell_inputs.counts_ch,
            ch_percell_inputs.solodir_ch,
            ch_percell_inputs.logs_ch,
            ch_percell_inputs.cutoff_ch,
            ch_percell_inputs.filtered_ch
        )
        percell_json = PERCELL_METRICS.out.percell_json

        // Build anchor from whichever mapping software was run
        def ch_anchor = star_summaries.mix(af_meta_info)

        // quant.json's num_genes counts USA matrix columns (3 per gene)
        def ch_af_mat_cols = af_mtx.map { meta, dir ->
            [ meta, file("${dir}/quants_mat_cols.txt") ]
        }

        // Join channels that need renaming by meta ID. af_meta_info is not joined again:
        // the anchor already carries it, and a second copy collides when staged.
        ch_to_rename = ch_anchor
            .join(cell_stats,     remainder: true)
            .join(knee_files,     remainder: true)
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
        def ch_analytical_manifest = ch_star_logs
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
            ch_star_logs.map{ it[1] }.collect().ifEmpty([]),
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
            antisense_metrics.collect().ifEmpty([]),
            secondderiv_knee.map { it[1] }.collect().ifEmpty([]),
            secondderiv_stats.map { it[1] }.collect().ifEmpty([]),
            percell_json.collect().ifEmpty([]),
            cs_ambient_hist_plot.collect().ifEmpty([]),
            cs_umap_comparison_plot.collect().ifEmpty([]),
            cs_top_genes.collect().ifEmpty([]),
            geneext_report.collect().ifEmpty([]),
            geneext_log.collect().ifEmpty([])
        )

        // Every file the mapping statistics scripts read, paired with the path it has
        // under the outdir. MAPPING_STATS rebuilds that layout from the staged files,
        // so it starts once they all exist rather than reading the outdir, which
        // publishDir fills asynchronously.
        //
        // The path helpers are invoked with .call(): Nextflow 26.04 does not resolve a closure
        // variable called like a function. ch_mapping_stats_files is assigned without 'def', as
        // Nextflow 24.04 fails to compile this declaration with it.
        def recalled  = params.cellfilter_method in ["second_derivative", "manual_cutoff"]
        def star_dir  = { meta -> "mapping_STARsolo/${meta.id}" }
        def gene_dir  = { meta -> "mapping_STARsolo/${meta.id}/${meta.id}_Solo.out/GeneFull_Ex50pAS" }
        def af_counts = { meta -> "mapping_alevin/${meta.id}/${meta.id}_counts" }

        ch_mapping_stats_files = ch_star_logs.map { meta, f -> ["${star_dir.call(meta)}/${f.name}", f] }
            .mix(
                star_full_logs.map  { meta, f -> ["${star_dir.call(meta)}/${f.name}", f] },
                star_summaries.map  { meta, f -> ["${gene_dir.call(meta)}/Summary.csv", f] },
                cell_stats.map      { meta, f -> ["${gene_dir.call(meta)}/CellReads.stats", f] },
                starsolo_raw_mtx.map { meta, dir -> ["${gene_dir.call(meta)}/raw", dir] },

                // The cell set: FILTER_MATRICES' matrix when cells are re-called, STARsolo's own otherwise
                star_filtered_mtx.map { meta, dir ->
                    [recalled ? "${gene_dir.call(meta)}/filtered_secondderiv/filtered" : "${gene_dir.call(meta)}/filtered", dir]
                },
                ch_star_logs.join(secondderiv_stats).map { meta, _log, f ->
                    ["${gene_dir.call(meta)}/filtered_secondderiv/${f.name}", f]
                },

                af_meta_info.map   { meta, f -> ["mapping_alevin/${meta.id}/${meta.id}_run/aux_info/meta_info.json", f] },
                af_quant_json.map  { meta, f -> ["${af_counts.call(meta)}/quant.json", f] },
                af_cell_meta.map   { meta, f -> ["${af_counts.call(meta)}/cell_meta.tsv", f] },
                ch_af_mat_cols.map { meta, f -> ["${af_counts.call(meta)}/alevin/quants_mat_cols.txt", f] },
                af_meta_info.join(secondderiv_stats).map { meta, _info, f ->
                    ["${af_counts.call(meta)}/alevin/filtered_secondderiv/${f.name}", f]
                },

                vendor_results.map { pipeline, meta, f ->
                    pipeline == 'cellranger'
                        ? ["CellRanger_pipeline/${meta.id}_count/outs", f]
                        : f.name == 'sample_all_stats.csv'
                            ? ["ParseBio_pipeline/${meta.id}/all-sample/report/${f.name}", f]
                            : ["ParseBio_pipeline/${meta.id}/${f.name}", f]
                }
            )
            // Sorted, so the staged inputs hash the same on every -resume
            .toSortedList { a, b -> a[0].toString() <=> b[0].toString() }
            .map { rows -> [ rows.collect { row -> row[0].toString() }, rows.collect { row -> row[1] } ] }

        MAPPING_STATS(ch_mapping_stats_files)

        // The outputs MultiQC has a module for: FastQC zips, STAR's Log.final.out
        // and the Kraken2 reports
        def ch_multiqc_files = fastqc_results
            .flatten()
            .filter { f -> f.name.endsWith('_fastqc.zip') }
            .mix(
                ch_star_logs.map { meta, log_final -> log_final },
                kraken_report
            )
            .collect()
            .ifEmpty([])

        // Salmon's meta_info.json takes its sample name from the run directory it sits
        // in, so the IDs are passed along to rebuild '<id>_run/aux_info/' in MULTIQC
        def ch_multiqc_salmon = af_meta_info
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .map { rows -> [ rows.collect { row -> row[0].id }, rows.collect { row -> row[1] } ] }

        ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        MULTIQC(
            ch_multiqc_files,
            ch_multiqc_salmon,
            ch_multiqc_config
        )

        def ch_versions = PERCELL_METRICS.out.versions.mix(
            PREPARE_DASHBOARD_INPUTS.out.versions,
            GENERATE_DASHBOARD.out.versions,
            MAPPING_STATS.out.versions,
            MULTIQC.out.versions
        )

    emit:
        dashboard_html  = GENERATE_DASHBOARD.out.html
        multiqc_report  = MULTIQC.out.report
        versions        = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
