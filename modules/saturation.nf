process SATURATION {
    publishDir "${params.resDir}/saturation/${sample_id}", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(mapping_files)
    file(bam_index)

    output:
    path("saturation*")

    script:
    """
    echo "\n\n==================  SATURATION =================="
    echo "Processing files: ${mapping_files}"

    # Find the correct files from the list (mapping_files)
    summary_file=\$(ls *Solo.out/Gene/Summary.csv | head -n 1)
    log_final_file=\$(ls *Log.final.out | head -n 1)
    bam_file=\$(ls ${sample_id}_Aligned.sortedByCoord.out.bam | head -n 1)

    echo "Summary file: \${summary_file}"
    echo "Log final file: \${log_final_file}"
    echo "BAM file: \${bam_file}"

    n_cells=\$( cat \${summary_file} | grep 'Estimated Number of Cells' | sed 's/,/ /g' | awk '{print \$NF}' )
    n_reads=\$( cat \${log_final_file} | grep 'Number of input reads' | awk '{print \$NF}' )
    MAPREADS=\$( samtools view -F 260 \${bam_file} | wc -l )
    map_rate=\$( echo "scale=4; \${MAPREADS}/\${n_reads}" | bc )
    temp_folder="_tmp_${sample_id}"
    echo "cells:\${n_cells} reads:\${n_reads} mapreads:\${MAPREADS} maprate:\${map_rate}"

    python ${params.baseDir}/submodules/10x_saturate/saturation_table.py \\
            --bam \${bam_file} \\
            --ncells \${n_cells} \\
            --mapping_rate \${map_rate} \\
            --temp \${temp_folder} \\
            --output saturation_output.tsv 

    python ${params.baseDir}/submodules/10x_saturate/scripts/plot_curve.py  \\
            saturation_output.tsv \\
            saturation.png \\
            --target 0.7 \\
            > saturation.log 2>&1   

    """
}
