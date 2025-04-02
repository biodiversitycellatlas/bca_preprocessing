process SATURATION_PLOT {
    publishDir "${params.resDir}/saturation/${config}/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), path(mapping_files)
    file(saturation_tsv)
    val(config)

    output:
    path("saturation*")

    script:
    """
    echo "\n\n==================  SATURATION PLOT =================="
    echo "Processing files: ${mapping_files}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    python ${params.baseDir}/submodules/10x_saturate/scripts/plot_curve.py  \\
            ${saturation_tsv} \\
            saturation.png \\
            --target 0.7 \\
            > saturation.log 2>&1   

    """
}
