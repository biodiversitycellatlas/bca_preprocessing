// ====================  SATURATION  ===================== \\ 
// Creates saturation plots using the tool: 10x_saturate   \\
// which is an external package linked using the github    \\
// submodule function.                                     \\

process SATURATION {
    publishDir "${params.resDir}/saturation/saturation_plots_${config_name}/${sample_id}", mode: 'symlink'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)
    file(bam_index)

    output:
    path("*")

    script:
    """
    echo "\n\n==================  SATURATION ${config_name} =================="
    # Verify the input files
    echo "Processing files: ${mapping_files}"

    # Find the correct files from the list (mapping_files)
    summary_file=\$(ls *Solo.out/Gene/Summary.csv | head -n 1)
    log_final_file=\$(ls *Log.final.out | head -n 1)
    bam_file=\$(ls *Aligned.sortedByCoord.out.bam | head -n 1)

    echo "Summary file: \${summary_file}"
    echo "Log final file: \${log_final_file}"
    echo "BAM file: \${bam_file}"

    n_cells=\$( cat \${summary_file} | grep 'Estimated Number of Cells' | sed 's/,/ /g' | awk '{print \$NF}' )
    n_reads=\$( cat \${log_final_file} | grep 'Number of input reads' | awk '{print \$NF}' )
    MAPREADS=\$( samtools view -F 260 \${bam_file} | wc -l )
    map_rate=\$( echo "scale=4; \${MAPREADS}/\${n_reads}" | bc )
    temp_folder="_tmp_${sample_id}_${config_name}"
    echo "cells:\${n_cells} reads:\${n_reads} mapreads:\${MAPREADS} maprate:\${map_rate}"

    python ${params.baseDir}/ext_programs/10x_saturate/saturation_table.py \
            --bam \${bam_file} \
            --ncells \${n_cells} \
            --mapping_rate \${map_rate} \
            --temp \${temp_folder} \
            --output output.tsv 
    # --code_dir "${params.baseDir}/ext_programs/10x_saturate/scripts"

    python ${params.baseDir}/ext_programs/10x_saturate/scripts/plot_curve.py  \
            output.tsv \
            saturation.png \
            --target 0.7    

    """
}
