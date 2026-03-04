process GENERATE_DASHBOARD {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label 'process_single'

    input:
    path(samplesheet)
    path(run_config)
    path(star_logs)
    path(star_summaries)
    path(star_full_logs)
    path(saturation_logs)
    path(cell_stats)
    path(af_meta_info)
    path(af_quant_json)
    path(af_cell_meta)
    path(sankey_files)
    path(saturation_imgs)
    path(residuals_imgs)
    path(knee_files)
    path(mt_rrna_metrics)
    path(per_cell_files)

    output:
    path "dashboard.html"  , emit: html
    path "versions.yml"    , emit: versions

    when:
    task.afterScript = 'echo "Dashboard generation complete"'

    script:
    def args = task.ext.args ?: ''
    """
    echo "======== Generating Dashboard ==========="
    echo "Samplesheet: ${samplesheet}"
    echo "Run config: ${run_config}"
    echo "STAR logs: ${star_logs}"
    echo "STAR summaries: ${star_summaries}"
    echo "STAR full logs: ${star_full_logs}"
    echo "Saturation logs: ${saturation_logs}"
    echo "Cell stats: ${cell_stats}"
    echo "Alevin-fry meta_info: ${af_meta_info}"
    echo "Alevin-fry quant.json: ${af_quant_json}"
    echo "Alevin-fry cell_meta: ${af_cell_meta}"
    echo "Sankey files: ${sankey_files}"
    echo "Saturation images: ${saturation_imgs}"
    echo "Residuals images: ${residuals_imgs}"
    echo "Knee files: ${knee_files}"
    echo "Mitochondrial rRNA metrics: ${mt_rrna_metrics}"
    echo "Per-cell metrics: ${per_cell_files}"

    # Execute the python dashboard generator
    generate_dashboard.py \\
        --project "bca_preprocessing" \\
        --samplesheet ${samplesheet} \\
        --template ${projectDir}/bin/dashboard_report.html \\
        --run_config ${run_config} \\
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
        --sankey_files ${sankey_files} \\
        --per_cell_files ${per_cell_files} \\
        --saturation_imgs ${saturation_imgs} \\
        --residuals_imgs ${residuals_imgs} \\
        --knee_files ${knee_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
