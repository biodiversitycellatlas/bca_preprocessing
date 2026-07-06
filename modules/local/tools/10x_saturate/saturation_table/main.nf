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

    output:
    tuple val(meta), path("saturation_output.tsv"), emit: saturation_table
    path "versions.yml",                            emit: versions

    script:
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
