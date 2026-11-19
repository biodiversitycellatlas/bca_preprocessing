process SCDBLFINDER {
    publishDir path: { "${params.outdir}/doublet_filtering/${meta.id}/scdblfinder" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${meta.mapping_method}"

    // Memory tracks the cell-called matrix
    cpus   1
    time   { 30.m * task.attempt }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.SCDBLFINDER, [tenx_dir], task.attempt, 10) }

    // scDblFinder cannot model every matrix: its kNN step needs more cells than a poor
    // library retains. A failure never ends the run (error_optional: one retry when the
    // task was killed for its resources, otherwise ignored); the sample's doublets then
    // come from Scrublet alone, and filtering_workflow.nf says so.
    label 'error_optional'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_bioconductor-scdblfinder_r-xgboost_r-argparse_pruned:60e56785b2c0e166' :
        'community.wave.seqera.io/library/r-base_bioconductor-scdblfinder_r-xgboost_r-argparse_pruned:1a632f336ca0f064' }"

    input:
    tuple val(meta), path(tenx_dir)

    output:
    tuple val(meta), path("${meta.id}_scDblFinder_results.tsv"), emit: scdblfinder_results
    path("${meta.id}_scDblFinder_summary.tsv"),                  emit: scdblfinder_summary, optional: true
    path "versions.yml",                                         emit: versions

    script:
    """
    echo -e "\\n\\n==================  scDblFinder =================="
    echo "Meta: ${meta}"
    echo "10x dir: ${tenx_dir}"
    echo "Datatype: ${meta.datatype}"

    demuxafy_scdblfinder.R \\
        -o scdblfinder_out \\
        -t ${tenx_dir} \\
        --doublet_rate_per_1000 ${params.doublet_rate_per_1000}

    mv scdblfinder_out/scDblFinder_doublets_singlets.tsv ${meta.id}_scDblFinder_results.tsv
    if [ -f scdblfinder_out/scDblFinder_doublet_summary.tsv ]; then
        mv scdblfinder_out/scDblFinder_doublet_summary.tsv ${meta.id}_scDblFinder_summary.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        scDblFinder: \$(Rscript -e 'cat(as.character(packageVersion("scDblFinder")))')
    END_VERSIONS
    """
}
