process PAVIAN {
    publishDir path: { "${params.outdir}/kraken" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_single2'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/r-base_r-dplyr_r-readr_r-purrr_pruned:612296ee176e80a7"
    // The script lives in a submodule. A container only sees the task directory and bin/,
    // unless the pipeline happens to share a parent directory with the work directory
    containerOptions { "${workflow.containerEngine in ['docker', 'podman'] ? '-v' : '-B'} ${projectDir}/submodules/pavianCore:${projectDir}/submodules/pavianCore" }

    input:
    path(kraken_out)

    output:
    path("*.sankey.html"), emit: sankey
    path "versions.yml",   emit: versions

    script:
    // sankeyD3 is only on GitHub, so neither the conda env nor the container ships it; its
    // dependencies do come with both. Under conda it is installed into the env the first time
    // PAVIAN runs, and the lock keeps parallel PAVIAN tasks from installing it at the same
    // time. A container image is read-only, so there it goes into the task directory instead,
    // where no lock is needed.
    """
    install_sankey='if (!requireNamespace("sankeyD3", quietly = TRUE)) remotes::install_github("fbreitwieser/sankeyD3", upgrade = "never", dependencies = NA)'
    if [ -n "\${CONDA_PREFIX:-}" ] && [ -w "\${CONDA_PREFIX}" ]; then
        ( flock 9; Rscript --vanilla -e "\$install_sankey" ) 9> "\${CONDA_PREFIX}/.sankeyD3.lock"
    else
        export R_LIBS_USER="\$PWD/R_libs"
        mkdir -p "\$R_LIBS_USER"
        Rscript --vanilla -e "\$install_sankey"
    fi

    Rscript ${projectDir}/submodules/pavianCore/exec/pavianCoreTools.R \\
        --input ${kraken_out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
