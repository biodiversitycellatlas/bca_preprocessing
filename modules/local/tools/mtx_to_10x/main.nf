process MTX_TO_10X {
    publishDir "${params.outdir}/doublet_filtering/${meta.id}/10x_export", mode: 'copy'
    tag "${meta.id} | ${mapping_method}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(mtx), path(barcodes), path(features), val(mapping_method), val(datatype)

    output:
    tuple val(meta), path("${meta.id}_10x"), val(mapping_method), val(datatype), emit: tenx_dir

    script:
    """
    echo "\n\n==================  MTX to 10x export =================="
    echo "Meta: ${meta}"
    echo "Mapping method: ${mapping_method}"

    mtx_to_10x.py \\
        --mtx ${mtx} \\
        --barcodes ${barcodes} \\
        --features ${features} \\
        --outdir ${meta.id}_10x
    """
}
