process SATURATION_TABLE {
    publishDir "${params.outdir}/saturation/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/pysam_samtools_bc_python_pruned:82a1e27e868113f0"

    input:
    tuple val(meta), file(bam_file)
    tuple val(meta), file(star_summary_file)
    tuple val(meta), file(star_log_final_file)
    tuple val(meta), file(samtools_bai)
    tuple val(meta), file(samtools_mapreads)
    tuple val(meta), file(secondderiv_stats)

    output:
    tuple val(meta), path("saturation_output.tsv"), emit: saturation_table
    path "versions.yml",                            emit: versions

    script:
    // Present only for cellfilter_method = "second_derivative"
    def sd_stats_file = secondderiv_stats ?: ''
    """
    echo "\n\n==================  SATURATION TABLE =================="
    echo "BAM file: ${bam_file}"
    echo "BAM index: ${samtools_bai}"
    echo "Mapped reads: ${samtools_mapreads}"
    echo "Summary file: ${star_summary_file}"
    echo "Log final file: ${star_log_final_file}"

    # Read the mapped reads from the file
    MAPREADS=\$( cat ${samtools_mapreads} )

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
    temp_folder="_tmp"

    echo "cells:\${n_cells} reads:\${n_reads} mapreads:\${MAPREADS} maprate:\${map_rate}"

    python ${projectDir}/submodules/10x_saturate/saturation_table.py \\
        --bam ${bam_file} \\
        --ncells \${n_cells} \\
        --mapping_rate \${map_rate} \\
        --temp \${temp_folder} \\
        --output saturation_output.tsv
    echo "Created saturation_output.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
