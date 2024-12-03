// ====================  SATURATION  ===================== \\ 
// Creates saturation plots using the tool: 10x_saturate   \\
// which is an external package linked using the github    \\
// submodule function.                                     \\

process SATURATION {
    publishDir "${params.resDir}/saturation_plots_${config_name}", mode: 'symlink'
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    output:
    path("*")

    script:
    """
    echo "\n\n==================  SATURATION ${config_name} =================="
    # Verify the input files
    echo "Processing bam file in Dir: ${mapping_files}"
    echo "Directory name: ${mapping_files.baseName}"

    n_cells=\$( cat ${mapping_files}/Solo.out/GeneFull/Summary.csv | grep 'Estimated Number of Cells' | sed 's/,/ /g' | awk '{print \$NF}' )
    n_reads=\$( cat ${mapping_files}/Log.final.out | grep 'Number of input reads' | awk '{print \$NF}' )
    MAPREADS=\$( samtools view -F 260 ${mapping_files}/*.bam | wc -l )
    map_rate=\$( echo "scale=4; \${MAPREADS}/\${n_reads}" | bc )
    temp_folder="${params.codeDir}/_tmp_${mapping_files.baseName}"
    echo "cells:\${n_cells} reads:\${n_reads} mapreads:\${MAPREADS} maprate:\${map_rate}"

    python ${params.codeDir}/ext_programs/10x_saturate/saturation_table.py \
            -b \${mapping_files}/*.bam \
            -n \${n_cells} \
            -r \${map_rate} \
            -t \${temp_folder} \
            -o output.tsv

    python ${params.codeDir}/ext_programs/10x_saturate/scripts/plot_curve.py  \
            output.tsv \
            saturation.png \
            --target 0.7

    """
}
