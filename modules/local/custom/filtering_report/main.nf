process GENERATE_FILTERING_REPORT {
    publishDir path: { "${params.outdir}" }, mode: 'copy', overwrite: true, saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_single_mem2'

    container "oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930"
    conda "${moduleDir}/environment.yml"

    input:
    path(summaries)
    path(mt_rrna_metrics)
    path(antisense_metrics)
    path(template)
    path(logo)
    val(expected_cells)
    val(thresholds)

    output:
    path "filtering_report.html", emit: html
    path "versions.yml",          emit: versions

    script:
    // A threshold left null starts switched off in the report
    def flags = [
        max_mito_pct    : '--max_mito_pct',
        max_doublet_pct : '--max_doublet_pct',
        max_alpha_hat   : '--max_alpha_hat',
        min_gene_umis   : '--min_gene_umis',
    ]
    def threshold_args = flags
        .findAll { key, _flag -> thresholds[key] != null }
        .collect { key, flag -> "${flag} ${thresholds[key]}" }
        .join(' ')
    """
    generate_filtering_report.py \\
        --template ${template} \\
        --summaries ${summaries} \\
        --mt_rrna_metrics ${mt_rrna_metrics} \\
        --antisense_metrics ${antisense_metrics} \\
        --expected_cells ${expected_cells.join(' ')} \\
        ${threshold_args} \\
        --pfam_presets '${thresholds.pfam_presets ?: ''}' \\
        --version "${workflow.manifest.version}" \\
        --logo ${logo} \\
        --output filtering_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch filtering_report.html versions.yml
    """
}
