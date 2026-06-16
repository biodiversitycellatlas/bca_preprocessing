process PAVIAN {
    publishDir "${params.outdir}/kraken2", mode: 'copy'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/r-base:4.5.2--d0821b0cfd27eb07"

    input:
    path(kraken_out)

    output:
    path("*.sankey.html")

    script:
    """
    # Install pavianCoreTools packages
    chmod +x ${projectDir}/submodules/pavianCore/exec/install_pavianCoreTools_packges.R
    Rscript ${projectDir}/submodules/pavianCore/exec/install_pavianCoreTools_packges.R

    # Make pavianCoreTools script executable
    chmod +x ${projectDir}/submodules/pavianCore/exec/pavianCoreTools.R

    # Run pavianCoreTools
    Rscript ${projectDir}/submodules/pavianCore/exec/pavianCoreTools.R \\
        --input ${kraken_out}
    """
}
