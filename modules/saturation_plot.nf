process SATURATION_PLOT {
    publishDir "${params.output_dir}/saturation/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), path(mapping_files)
    file(saturation_tsv)

    output:
    path("saturation*")

    script:
    """
    echo "\n\n==================  SATURATION PLOT =================="
    echo "Processing files: ${mapping_files}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    python ${params.code_dir}/submodules/10x_saturate/scripts/plot_curve.py  \\
            ${saturation_tsv} \\
            saturation.png \\
            --target 0.7 \\
            > saturation.log 2>&1   

    """
}
