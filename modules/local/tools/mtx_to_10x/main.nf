process MTX_TO_10X {
    publishDir "${params.outdir}/doublet_filtering/${meta.id}/10x_export", mode: 'copy'
    tag "${meta.id} | ${meta.mapping_method}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(mtx), path(barcodes), path(features)

    output:
    tuple val(meta), path("${meta.id}_10x"), emit: tenx_dir

    script:
    """
    echo "\n\n==================  MTX to 10x export =================="
    echo "Meta: ${meta}"

    mtx_to_10x.py \\
        --mtx ${mtx} \\
        --barcodes ${barcodes} \\
        --features ${features} \\
        --outdir ${meta.id}_10x
    """
}
