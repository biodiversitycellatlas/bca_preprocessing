// ==================  MAPPING STATISTICS  =================== \\ 

process MAPPING_STATS {
    publishDir "${params.resDir}/mapping_stats", mode: 'copy'
    debug true

    input:
    tuple val(sample_id_N), val(config_name_N), path(files_N)
    tuple val(sample_id_CR), val(config_name_CR), path(files_CR)

    output:
    file("star_mapping_stats.tsv")

    script:
    """
    echo "\n\n===============  MAPPING STATISTICS  ==============="

    # Change to the directory that contains the "data" directory
    cd ${params.resDir} || { echo "Error: Could not cd to script directory"; exit 1; }

    # Define output file
    output_file="star_mapping_stats.tsv"

    # Print header to output file
    echo -e "Directory\\tSample\\tN reads/sample\\tN R1 >Q30\\tN R2 >Q30\\tN uniquely mapped reads\\t% uniquely mapped reads\\t% multi-mapped reads\\t% multi-mapped reads: too many\\t% unmapped: too short\\t% unmapped: other\\tN cells\\tUMI cutoff used for cell calling\\tsaturation\\tNoise (% UMIs in non-cell barcodes)\\tMedian % rRNA\\tMedian % mtDNA" > "\$output_file"

    # For each mapping directory under the experiment
    for map_dir in ${params.resDir}/mapping_STARsolo_*; do
        [ -d "\$map_dir" ] || continue
        # For each sample directory under the mapping directory
        for sample_dir in "\$map_dir"/*; do
            [ -d "\$sample_dir" ] || continue

            LOG="\$sample_dir/Log.final.out"
            [ -f "\$LOG" ] || continue

            # Extract from STAR Log
            n_reads=\$(grep "Number of input reads" "\$LOG" | awk '{print \$NF}')
            n_uniquely_mapped=\$(grep "Uniquely mapped reads number" "\$LOG" | awk '{print \$NF}')
            p_uniquely_mapped=\$(grep "Uniquely mapped reads %" "\$LOG" | awk '{print \$NF}')
            p_multi_mapped=\$(grep "% of reads mapped to multiple loci" "\$LOG" | awk '{print \$NF}')
            p_multi_too_many=\$(grep "% of reads mapped to too many loci" "\$LOG" | awk '{print \$NF}')
            p_unmapped_short=\$(grep "% of reads unmapped: too short" "\$LOG" | awk '{print \$NF}')
            p_unmapped_other=\$(grep "% of reads unmapped: other" "\$LOG" | awk '{print \$NF}')

            # Extract QC metrics for R1 and R2 > Q30 (adjust as needed)
            # Here we assume these files contain just the numeric value.
            # If not found, default to "NA"
            R1_Q30=\$(cat "\$sample_dir/R1_q30.txt" 2>/dev/null || echo "NA")
            R2_Q30=\$(cat "\$sample_dir/R2_q30.txt" 2>/dev/null || echo "NA")

            # Check Gene and GeneFull directories in Solo.out
            solo_gene_dir="\$sample_dir/Solo.out/Gene/filtered"
            solo_genefull_dir="\$sample_dir/Solo.out/GeneFull/filtered"

            if [ -d "\$solo_gene_dir" ]; then
                solo_dir="\$solo_gene_dir"
            elif [ -d "\$solo_genefull_dir" ]; then
                solo_dir="\$solo_genefull_dir"
            else
                solo_dir=""
            fi

            # N cells (if barcodes.tsv is available)
            if [ -n "\$solo_dir" ] && [ -f "\$solo_dir/barcodes.tsv" ]; then
                N_cells=\$(wc -l < "\$solo_dir/barcodes.tsv")
            else
                N_cells="NA"
            fi

            # Extract other metrics from summary file (if available)
            if [ -n "\$solo_dir" ] && [ -f "\$solo_dir/summary.tsv" ]; then
                UMI_cutoff=\$(grep "UMI_cutoff" "\$solo_dir/summary.tsv" | cut -f2 || echo "NA")
                saturation=\$(grep "Saturation" "\$solo_dir/summary.tsv" | cut -f2 || echo "NA")
                noise=\$(grep "Noise_fraction" "\$solo_dir/summary.tsv" | cut -f2 || echo "NA")
                median_rRNA=\$(grep "Median_percent_rRNA" "\$solo_dir/summary.tsv" | cut -f2 || echo "NA")
                median_mtDNA=\$(grep "Median_percent_mtDNA" "\$solo_dir/summary.tsv" | cut -f2 || echo "NA")
            else
                UMI_cutoff="NA"
                saturation="NA"
                noise="NA"
                median_rRNA="NA"
                median_mtDNA="NA"
            fi

            # Append results to the output file
            echo -e "${params.resDir}\\t\$(basename "\$sample_dir")\\t\$n_reads\\t\$R1_Q30\\t\$R2_Q30\\t\$n_uniquely_mapped\\t\$p_uniquely_mapped\\t\$p_multi_mapped\\t\$p_multi_too_many\\t\$p_unmapped_short\\t\$p_unmapped_other\\t\$N_cells\\t\$UMI_cutoff\\t\$saturation\\t\$noise\\t\$median_rRNA\\t\$median_mtDNA" >> "\$output_file"
        done
    done


    echo "Finished. Results are in \$output_file"

    """
}
