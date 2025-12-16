process SATURATION_TABLE {
    publishDir "${params.outdir}/saturation/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/pysam_samtools_matplotlib_numpy_pruned:b8f551e4a5153343"

    input:
    tuple val(meta), path(mapping_files)
    file(filtered_bam)
    file(bam_index)
    val MAPREADS

    output:
    path("saturation_output.tsv")

    script:
    """
    echo "\n\n==================  SATURATION TABLE =================="
    echo "Processing files: ${mapping_files}"
    echo "Filtered BAM file: ${filtered_bam}"
    echo "Filtered BAM index file: ${bam_index}"
    echo "Mapped reads: ${MAPREADS}"

    # Find the correct files from the list (mapping_files)
    summary_file=\$(ls *Solo.out/Gene/Summary.csv | head -n 1)
    log_final_file=\$(ls *Log.final.out | head -n 1)
    bam_file=\$(ls ${meta.id}_Aligned.filtered.sorted.bam | head -n 1)

    echo "Summary file: \${summary_file}"
    echo "Log final file: \${log_final_file}"
    echo "BAM file: \${bam_file}"

    n_cells=\$( cat \${summary_file} | grep 'Estimated Number of Cells' | sed 's/,/ /g' | awk '{print \$NF}' )
    n_reads=\$( cat \${log_final_file} | grep 'Number of input reads' | awk '{print \$NF}' )

    map_rate=\$( echo "scale=4; ${MAPREADS}/\${n_reads}" | bc )
    temp_folder="_tmp"
    echo "cells:\${n_cells} reads:\${n_reads} mapreads:${MAPREADS} maprate:\${map_rate}"

    python ${projectDir}/submodules/10x_saturate/saturation_table.py \\
        --bam \${bam_file} \\
        --ncells \${n_cells} \\
        --mapping_rate \${map_rate} \\
        --temp \${temp_folder} \\
        --output saturation_output.tsv
    echo "Created saturation_output.tsv"
    """
}
