process MAPPING_STATS {
    publishDir path: { "${params.outdir}/summary_results" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_single_mem2'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/pandas_python_r-base_r-data.table_pruned:30aa8e8a69298cb8"

    input:
    // Each staged file with the path it has under the outdir, in the same order
    tuple val(dests), path(files, stageAs: "staged_inputs/?/*")

    output:
    path("mapping_stats.tsv"),  emit: stats
    path("R_images"),           emit: plots, optional: true
    path "versions.yml",        emit: versions

    script:
    // Named after the outdir, since plot_umidist_cellgenecount.R labels its plots with it
    def root      = file(params.outdir).name ?: 'results'
    def file_list = files instanceof List ? files : [files]

    // Rebuild the outdir layout the scripts walk, e.g. mapping_STARsolo/<id>/<id>_Log.final.out
    def links = [dests, file_list].transpose().collect { dest, f ->
        def parent = dest.contains('/') ? dest.substring(0, dest.lastIndexOf('/')) : '.'
        "mkdir -p '${root}/${parent}' && ln -s \"\$PWD/${f}\" '${root}/${dest}'"
    }.join('\n    ')

    """
    mkdir -p '${root}'
    ${links}

    # sci-rocket results are not produced by this pipeline, so a run placed in the outdir is used as is
    if [ -d "${file(params.outdir)}/sci-rocket" ]; then
        ln -s "${file(params.outdir)}/sci-rocket" '${root}/sci-rocket'
    fi

    # Produce summary (.tsv) of mapping statistics
    dashboard_mappingstats.py '${root}'

    # UMI distribution and Cell + Gene count plots, for every run that mapped with STARsolo
    if [ -d '${root}/mapping_STARsolo' ]; then
        plot_umidist_cellgenecount.R '${root}'
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        r-base: \$(Rscript --version 2>&1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
