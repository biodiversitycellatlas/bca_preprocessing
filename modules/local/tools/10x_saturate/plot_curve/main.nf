process SATURATION_PLOT {
    publishDir path: { "${params.outdir}/saturation/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"
    label 'process_single2'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_samtools_matplotlib_numpy_pruned:b8f551e4a5153343' :
        'community.wave.seqera.io/library/pysam_samtools_matplotlib_numpy_pruned:c50e3d2d3c1c69f3' }"
    // The script lives in a submodule. A container only sees the task directory and bin/,
    // unless the pipeline happens to share a parent directory with the work directory
    containerOptions { "${workflow.containerEngine in ['docker', 'podman'] ? '-v' : '-B'} ${projectDir}/submodules/10x_saturate:${projectDir}/submodules/10x_saturate" }

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
    echo -e "\\n\\n==================  SATURATION PLOT =================="
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
