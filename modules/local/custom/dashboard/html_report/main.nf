process GENERATE_DASHBOARD {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label 'process_single'

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
    path(secondderiv_knee)
    path(secondderiv_stats)
    path(per_cell_files)
    path(cellsweep_plots_contrib)
    path(cellsweep_plots_umap)
    path(cellsweep_top_genes)

    output:
    path "dashboard.html"  , emit: html
    path "versions.yml"    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Execute the python dashboard generator
    generate_dashboard.py \\
        --project "bca_preprocessing" \\
        --samplesheet ${samplesheet} \\
        --template ${projectDir}/bin/dashboard_report.html \\
        --run_config ${run_config} \\
        --analytical_samples ${analytical_manifest} \\
        --output "dashboard.html" \\
        --star_logs ${star_logs} \\
        --star_summaries ${star_summaries} \\
        --star_full_logs ${star_full_logs} \\
        --mt_rrna_metrics ${mt_rrna_metrics} \\
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
        --cellsweep_tables ${cellsweep_top_genes} \\
        --cellsweep_plots_contrib ${cellsweep_plots_contrib} \\
        --cellsweep_plots_umap ${cellsweep_plots_umap}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
