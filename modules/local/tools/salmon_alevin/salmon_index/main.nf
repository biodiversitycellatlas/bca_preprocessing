process SALMON_INDEX {
    publishDir "${params.outdir}/genome", mode: 'copy'
    tag "${meta.id}"
    label 'process_high'


    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/salmon:1.10.3--726401738a281398"

    input:
    tuple val(meta), path(splici_index_reference)

    output:
    tuple val(meta), path("salmon_index_${meta.id}/"), emit: salmon_index


    script:
    """
    echo "\n\n==================  SALMON INDEX =================="
    # Build reference index
    salmon index \\
        -t ${splici_index_reference}/*.fa \\
        -i ./salmon_index_${meta.id} \\
        -k 31
    """
}
