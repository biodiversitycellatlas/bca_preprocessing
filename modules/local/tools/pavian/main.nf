process PAVIAN {
    publishDir "${params.outdir}/kraken", mode: 'copy'
    label 'process_single'

    conda "${projectDir}/work/envs/pavian"

    input:
    path(kraken_out)

    output:
    path("*.sankey.html")

    script:
    """
    Rscript ${projectDir}/submodules/pavianCore/exec/pavianCoreTools.R \\
        --input ${kraken_out}
    """
}
