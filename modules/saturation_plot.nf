process SATURATION_PLOT {
    publishDir "${params.output_dir}/saturation/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(mapping_files)
    file(saturation_tsv)

    output:
    path("saturation*")

    script:
    """
    echo "\n\n==================  SATURATION PLOT =================="
    echo "Processing files: ${mapping_files}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    python ${launchDir}/submodules/10x_saturate/scripts/plot_curve.py  \\
            ${saturation_tsv} \\
            saturation.png \\
            --target 0.7 \\
            > saturation.log 2>&1   

    """
}
