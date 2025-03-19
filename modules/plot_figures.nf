process PLOT_FIGURES {
    script:
    """
    echo "\n\n==================  PLOT FIGURES  =================="

    work_dir="${params.resDir}/mapping_STARsolo"
    
    # Create a CSV file with a header
    output_csv="\${work_dir}/intron_vs_mtdna_metrics.csv"
    echo "sample_name,config,mtDNA_fraction,fraction_intronic" > "\$output_csv"

    # Loop over subdirectories in $work_dir
    for config_dir in "\$work_dir"/*; do
        [[ ! -d "\$config_dir" ]] && continue  # Skip if not a directory
        echo "---- Plot Knee plots and Comparisons"
        python ${params.baseDir}/scripts/visualize_starsolo_comparison.py \${config_dir}

        echo "---- Plot Overview rRNA & mtDNA"
        python ${params.baseDir}/scripts/visualize_mt_rrna.py \${config_dir}

        # Loop over sample directories in each config directory: e.g. DPS in BCA001_..._DSP
        for sample_dir in "\${config_dir}"/*; do
            [[ ! -d "\${sample_dir}" ]] && continue
            echo "---- Percentage Mitochondrial DNA in Unique Reads"
            final_log="\${sample_dir}/Log.final.out"
            n_uniquely_mapped=\$(grep "Uniquely mapped reads number" "\${final_log}" | awk '{print \$NF}')

            sample_name=\$(basename \${sample_dir})
            config="\${config_dir##*_}"
            featcounts_dir="\${resDir}/rRNA_mtDNA/rRNA_mtDNA_\${config}/\${sample_name}"

            mtDNA_summary="\${featcounts_dir}/feat_counts_mtDNA.txt.summary"
            mtDNA_assigned=\$(grep "^Assigned" "\${mtDNA_summary}" | awk '{print \$2}')
            mtDNA_percentage=\$(awk -v sum=\${n_uniquely_mapped} -v assigned2=\${mtDNA_assigned} 'BEGIN{if(assigned2>0){printf("%.6f", (assigned2/sum))}else{print("NA")}}')

            echo "---- Percentage Intronic Reads"
            exonic_sum=\$(cat \${sample_dir}/Solo.out/Gene/raw/matrix.mtx | awk 'NR>3 {sum+=\$3} END{print sum}')
            fullgene_sum=\$(cat \${sample_dir}/Solo.out/GeneFull/raw/matrix.mtx | awk 'NR>3 {sum+=\$3} END{print sum}')
            intronic_sum=\$(( fullgene_sum - exonic_sum ))

            fraction_intronic=\$(awk -v i="\${intronic_sum}" -v t="\${fullgene_sum}" \
                                'BEGIN{printf("%.4f", (i)/t)}')

            # Append Metrics to CSV
            echo "\${sample_name},\${config},\${mtDNA_percentage},\${fraction_intronic}" >> "\$output_csv"
        done
    done

    echo "---- Plot Intronic Reads vs. Mitochondrial DNA %"
    python ${params.baseDir}/scripts/visualize_intronic_vs_mtdna.py "\$output_csv"

    """
}