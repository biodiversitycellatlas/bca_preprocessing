process SATURATION_PLOT {
    publishDir "${params.outdir}/saturation/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/pysam_samtools_matplotlib_numpy_pruned:b8f551e4a5153343"

    input:
    tuple val(meta), path(mapping_files)
    file(saturation_output)

    output:
    path("saturation*")

    script:
    """
    echo "\n\n==================  SATURATION PLOT =================="
    python ${projectDir}/submodules/10x_saturate/scripts/plot_curve.py  \\
        ${saturation_output} \\
        saturation.png \\
        --target 0.7 \\
        > saturation.log 2>&1 || true
    """
}
