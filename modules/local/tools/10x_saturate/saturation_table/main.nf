process SATURATION_TABLE {
    publishDir path: { "${params.outdir}/saturation/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"
    label 'error_retry'

    // Memory and time track the size of the BAM file
    cpus   4
    time   { BcaResources.scaledTime(
        params.dynamic_time?.SATURATION_TABLE, [bam_file], task.attempt, 4) }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.SATURATION_TABLE, [bam_file], task.attempt, 10) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_samtools_bc_python_pruned:82a1e27e868113f0' :
        'community.wave.seqera.io/library/pysam_samtools_bc_python_pruned:d34a762142fc1075' }"
        
    // The script lives in a submodule. A container only sees the task directory and bin/,
    // unless the pipeline happens to share a parent directory with the work directory
    containerOptions { "${workflow.containerEngine in ['docker', 'podman'] ? '-v' : '-B'} ${projectDir}/submodules/10x_saturate:${projectDir}/submodules/10x_saturate" }

    input:
    tuple val(meta),   path(bam_file)
    tuple val(_meta2), path(star_summary_file)
    tuple val(_meta3), path(star_log_final_file)
    tuple val(_meta4), path(bam_index)
    tuple val(_meta5), path(secondderiv_stats)

    output:
    tuple val(meta), path("saturation_output.tsv"), emit: saturation_table
    path "versions.yml",                            emit: versions

    script:
    // Present only for cellfilter_method = "second_derivative"
    def sd_stats_file = secondderiv_stats ?: ''
    """
    echo -e "\\n\\n==================  SATURATION TABLE =================="
    echo "BAM file: ${bam_file}"
    echo "BAM index: ${bam_index}"
    echo "Summary file: ${star_summary_file}"
    echo "Log final file: ${star_log_final_file}"

    # Mapped reads, counted once per read: uniquely mapped plus multimapped. This is the
    # number of primary mapped records in the BAM, which gettags.py keeps for the tags.
    MAPREADS=\$( awk -F'|' '/Uniquely mapped reads number|Number of reads mapped to multiple loci/ { gsub(/[ \\t]/, "", \$2); s += \$2 } END { print s + 0 }' ${star_log_final_file} )
    if [ "\${MAPREADS}" -eq 0 ]; then
        echo "[ERROR] no mapped reads in ${star_log_final_file}; cannot compute the mapping rate" >&2
        exit 1
    fi

    echo "Mapped reads: \${MAPREADS}"

    n_cells=\$( cat ${star_summary_file} | grep 'Estimated Number of Cells' | sed 's/,/ /g' | awk '{print \$NF}' )
    n_reads=\$( cat ${star_log_final_file} | grep 'Number of input reads' | awk '{print \$NF}' )

    # The saturation curve is fitted per cell, so it has to use the cell set the
    # matrices were filtered on. With cellfilter_method = "second_derivative" that is
    # the count recomputed by FILTER_MATRICES, not STARsolo's own Summary.csv estimate.
    if [ -s "${sd_stats_file}" ]; then
        sd_cells=\$(python3 -c "import json; print(json.load(open('${sd_stats_file}')).get('estimated_cells', ''))")
        if [ -n "\${sd_cells}" ]; then
            echo "Second-derivative cell count: \${sd_cells} (STARsolo estimated \${n_cells})"
            n_cells=\${sd_cells}
        else
            echo "[WARNING] no estimated_cells in ${sd_stats_file}; keeping STARsolo's estimate"
        fi
    fi

    map_rate=\$( echo "scale=4; \${MAPREADS}/\${n_reads}" | bc )

    # data.table otherwise sizes its thread pool from the node's cores, not the task's
    export R_DATATABLE_NUM_THREADS=${task.cpus}

    echo "cells:\${n_cells} reads:\${n_reads} mapreads:\${MAPREADS} maprate:\${map_rate}"

    # The steps of 10x_saturate's saturation_table.py, run here so each gets its own
    # thread count: that wrapper hands one --ncpu to both (default 1, so the tag dump ran
    # on one core) and ignores their exit codes. The tag dump runs samtools with -@ and
    # one worker per chromosome. The curve stays on one R process: mclapply would give
    # every fork its own subsample of the tag table and multiply the memory; data.table
    # still threads within it.
    sat_dir=${projectDir}/submodules/10x_saturate/scripts
    mkdir _tmp
    python \${sat_dir}/gettags.py ${bam_file} _tmp --num_threads ${task.cpus}
    cat _tmp/*.tsv > _tmp/tags.tab
    echo "Counting tags done: _tmp/tags.tab"

    Rscript \${sat_dir}/compute_saturation.r \\
        --ncpu 1 \\
        --mapping_rate \${map_rate} \\
        --nstep 20 \\
        _tmp/tags.tab \\
        saturation_output.tsv \\
        \${n_cells}
    echo "Created saturation_output.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
