process PAVIAN {
    publishDir "${params.outdir}/kraken", mode: 'copy'
    label 'process_single2'

    conda "${moduleDir}/environment.yml"

    input:
    path(kraken_out)

    output:
    path("*.sankey.html"), emit: sankey
    path "versions.yml",   emit: versions

    script:
    // sankeyD3 is only on GitHub, so it is installed into the conda env the first time
    // PAVIAN runs. The lock keeps parallel PAVIAN tasks from installing it at the same time.
    """
    (
        flock 9
        Rscript --vanilla -e '
        if (!requireNamespace("sankeyD3", quietly = TRUE)) {
            remotes::install_github("fbreitwieser/sankeyD3", upgrade = "never", dependencies = TRUE)
        }
        '
    ) 9> "\${CONDA_PREFIX:-.}/.sankeyD3.lock"

    Rscript ${projectDir}/submodules/pavianCore/exec/pavianCoreTools.R \\
        --input ${kraken_out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
