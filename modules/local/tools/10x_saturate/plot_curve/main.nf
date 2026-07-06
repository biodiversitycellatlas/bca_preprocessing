process SATURATION_PLOT {
    publishDir "${params.outdir}/saturation/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/pysam_samtools_matplotlib_numpy_pruned:b8f551e4a5153343"

    input:
    tuple val(meta), file(saturation_output)

    output:
    path("${meta.id}_saturation*"),                         emit: all
    path("${meta.id}_saturation.log"),                      emit: logs
    path("${meta.id}_saturation.png"),                      emit: img_saturation
    path("${meta.id}_saturation_residuals.png"),            emit: img_residuals
    path "versions.yml",                                    emit: versions

    script:
    """
    echo "\n\n==================  SATURATION PLOT =================="
    python ${projectDir}/submodules/10x_saturate/scripts/plot_curve.py  \\
        ${saturation_output} \\
        ${meta.id}_saturation.png \\
        --target ${params.saturation_target} \\
        > ${meta.id}_saturation.log 2>&1 || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
