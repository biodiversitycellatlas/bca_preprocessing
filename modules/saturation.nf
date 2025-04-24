process SATURATION {
    publishDir "${params.output_dir}/saturation/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(mapping_files)
    file(bam_index)

    output:
    file("saturation_output.tsv")

    script:
    """
    echo "\n\n==================  SATURATION =================="
    echo "Processing files: ${mapping_files}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    # Remove unmapped reads from the BAM file
    samtools view -b -F 4 ${meta.id}_Aligned.sortedByCoord.out.bam > ${meta.id}_Aligned.sortedByCoord.out.mapped.bam
    samtools index ${meta.id}_Aligned.sortedByCoord.out.mapped.bam

    # Find the correct files from the list (mapping_files)
    summary_file=\$(ls *Solo.out/Gene/Summary.csv | head -n 1)
    log_final_file=\$(ls *Log.final.out | head -n 1)
    bam_file=\$(ls ${meta.id}_Aligned.sortedByCoord.out.mapped.bam | head -n 1)

    echo "Summary file: \${summary_file}"
    echo "Log final file: \${log_final_file}"
    echo "BAM file: \${bam_file}"

    n_cells=\$( cat \${summary_file} | grep 'Estimated Number of Cells' | sed 's/,/ /g' | awk '{print \$NF}' )
    n_reads=\$( cat \${log_final_file} | grep 'Number of input reads' | awk '{print \$NF}' )
    MAPREADS=\$( samtools view -F 260 \${bam_file} | wc -l )
    map_rate=\$( echo "scale=4; \${MAPREADS}/\${n_reads}" | bc )
    temp_folder="_tmp"
    echo "cells:\${n_cells} reads:\${n_reads} mapreads:\${MAPREADS} maprate:\${map_rate}"

    python ${launchDir}/submodules/10x_saturate/saturation_table.py \\
            --bam \${bam_file} \\
            --ncells \${n_cells} \\
            --mapping_rate \${map_rate} \\
            --temp \${temp_folder} \\
            --output saturation_output.tsv 

    ls -lah
    echo "Created saturation_output.tsv"
    """
}
