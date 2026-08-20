process PAVIAN {
    publishDir "${params.outdir}/kraken", mode: 'copy'
    label 'process_single2'

    conda "${projectDir}/work/envs/pavian"

    input:
    path(kraken_out)

    output:
    path("*.sankey.html"), emit: sankey
    path "versions.yml",   emit: versions

    script:
    """
    Rscript ${projectDir}/submodules/pavianCore/exec/pavianCoreTools.R \\
        --input ${kraken_out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
