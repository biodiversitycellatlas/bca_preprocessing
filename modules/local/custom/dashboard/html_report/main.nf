process GENERATE_DASHBOARD {
    publishDir path: { "${params.outdir}" }, mode: 'copy', overwrite: true, saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    
    // The dashboard reads and embeds everything it is given, so memory tracks the total
    // size of its inputs, which grows with the number of runs. Coefficients live in
    // params.dynamic_memory; without an entry the request is process_single_mem2's 10 GB.
    cpus   1
    time   { 30.m * task.attempt }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.GENERATE_DASHBOARD,
        [samplesheet, run_config, analytical_manifest, star_logs, star_summaries, star_full_logs,
         saturation_logs, cell_stats, af_meta_info, af_quant_json, af_cell_meta, af_mat_cols,
         sankey_files, saturation_imgs, residuals_imgs, knee_files, mt_rrna_metrics,
         antisense_metrics, secondderiv_knee, secondderiv_stats, matrix_stats, per_cell_files,
         cellsweep_plots_contrib, cellsweep_plots_umap, cellsweep_top_genes, geneext_report,
         geneext_log],
        task.attempt, 10) }

    input:
    path(samplesheet)
    path(run_config)
    path(analytical_manifest)
    path(star_logs)
    path(star_summaries)
    path(star_full_logs)
    path(saturation_logs)
    path(cell_stats)
    path(af_meta_info)
    path(af_quant_json)
    path(af_cell_meta)
    path(af_mat_cols)
    path(sankey_files)
    path(saturation_imgs)
    path(residuals_imgs)
    path(knee_files)
    path(mt_rrna_metrics)
    path(antisense_metrics)
    path(secondderiv_knee)
    path(secondderiv_stats)
    path(matrix_stats)
    path(per_cell_files)
    path(cellsweep_plots_contrib)
    path(cellsweep_plots_umap)
    path(cellsweep_top_genes)
    path(geneext_report)
    path(geneext_log)

    output:
    path "dashboard.html"  , emit: html
    path "versions.yml"    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Execute the python dashboard generator
    generate_dashboard.py \\
        --project "bca_preprocessing" \\
        --version "${workflow.manifest.version}" \\
        --commit "${workflow.commitId ?: 'N/A'}" \\
        --samplesheet ${samplesheet} \\
        --template ${projectDir}/bin/dashboard_report.html \\
        --run_config ${run_config} \\
        --analytical_samples ${analytical_manifest} \\
        --output "dashboard.html" \\
        --star_logs ${star_logs} \\
        --star_summaries ${star_summaries} \\
        --star_full_logs ${star_full_logs} \\
        --mt_rrna_metrics ${mt_rrna_metrics} \\
        --antisense_metrics ${antisense_metrics} \\
        --saturation_logs ${saturation_logs} \\
        --cell_stats ${cell_stats} \\
        --af_meta_info ${af_meta_info} \\
        --af_quant_json ${af_quant_json} \\
        --af_cell_meta ${af_cell_meta} \\
        --af_mat_cols ${af_mat_cols} \\
        --sankey_files ${sankey_files} \\
        --per_cell_files ${per_cell_files} \\
        --saturation_imgs ${saturation_imgs} \\
        --residuals_imgs ${residuals_imgs} \\
        --knee_files ${knee_files} \\
        --secondderiv_knee ${secondderiv_knee} \\
        --secondderiv_stats ${secondderiv_stats} \\
        --matrix_stats ${matrix_stats} \\
        --cellsweep_tables ${cellsweep_top_genes} \\
        --cellsweep_plots_contrib ${cellsweep_plots_contrib} \\
        --cellsweep_plots_umap ${cellsweep_plots_umap} \\
        --geneext_report ${geneext_report} \\
        --geneext_log ${geneext_log}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
