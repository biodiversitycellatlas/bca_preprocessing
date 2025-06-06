#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=../logs/%x.%j.out
#SBATCH --error=../logs/%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --qos=vshort
#SBATCH --mem=6G
#SBATCH --job-name mapping_stats

# ------------------------------------------------------------------
# Define Variables
# ------------------------------------------------------------------
output_dir=$1
output_file="mapping_stats.tsv"

# Change to the directory that contains the data
cd ${output_dir} || { echo "Error: Could not cd to data directory"; exit 1; }

# Create the output file and print the header
echo -e "Software\tSample\tN reads/sample\tN R1 >Q30\tN R2 >Q30\tN uniquely mapped reads\t\
    % uniquely mapped reads\t% multi-mapped reads\t% multi-mapped reads: too many\t\
    % unmapped: too short\t% unmapped: other\tExpected % Doublets\tTarget N cells\t\
    N cells\tUMI cutoff used for cell calling\tSaturation\tReads for 0.7 saturation\t\
    Noise (% UMIs in non-cell barcodes)\t% Intronic reads\t% rRNA in Unique reads\t\
    % mtDNA in Unique reads\t%mtDNA in multimappers all pos\t%mtDNA in multimappers primary pos\t\
    Mean Reads per Cell\tMedian UMI Counts per Cell\tMedian Genes per Cell\tTotal Genes Detected" > "$output_file"


# ------------------------------------------------------------------
# For each STARsolo mapping directory
# ------------------------------------------------------------------

for map_dir in "${output_dir}/mapping_STARsolo"/*; do
    echo "mapping dir: ${map_dir}"

    # Extract sample name
    sample_name=$(basename ${map_dir})
    echo "sample: $sample_name"

    # Check for Log.final.out file and skip if not found
    LOG="$map_dir/${sample_name}_Log.final.out"

    # Extract metrics: Log.final.out file
    n_reads=$(grep "Number of input reads" "$LOG" | awk '{print $NF}')
    n_uniquely_mapped=$(grep "Uniquely mapped reads number" "$LOG" | awk '{print $NF}')
    p_uniquely_mapped=$(grep "Uniquely mapped reads %" "$LOG" | awk '{print $NF}')
    p_multi_mapped=$(grep "% of reads mapped to multiple loci" "$LOG" | awk '{print $NF}')
    p_multi_too_many=$(grep "% of reads mapped to too many loci" "$LOG" | awk '{print $NF}')
    p_unmapped_short=$(grep "% of reads unmapped: too short" "$LOG" | awk '{print $NF}')
    p_unmapped_other=$(grep "% of reads unmapped: other" "$LOG" | awk '{print $NF}')

    # Check Gene and GeneFull_Ex50pAS directories in Solo.out
    solo_gene_dir="$map_dir/${sample_name}_Solo.out/Gene"
    solo_genefull_dir="$map_dir/${sample_name}_Solo.out/GeneFull_Ex50pAS"

    if [ -d "$solo_genefull_dir" ]; then
        solo_dir="$solo_genefull_dir"
        config="GeneFull_Ex50pAS"
    elif [ -d "$solo_gene_dir" ]; then
        solo_dir="$solo_gene_dir"
    fi
    echo "solo dir: $solo_dir"

    # Extract metrics: barcodes.csv file 
    N_cells=$(wc -l < "$solo_dir/filtered/barcodes.tsv")

    # Extract metrics: summary.csv file
    R1_Q30=$(grep "RNA" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
    R2_Q30=$(grep "CB+UMI" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
    saturation=$(grep "Saturation" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
    mean_reads=$(grep "Mean Reads per Cell" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
    median_umis=$(grep "Median UMI per Cell" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
    median_genes=$(grep "Median $config per Cell" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
    total_genes=$(grep "Total $config Detected" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")

    # Extract metrics: STAR Log file
    UMI_cutoff=$(grep "nUMImin" "$map_dir/${sample_name}_Log.out" | awk -F ';' '{print $2}' | head -n 1 | awk -F '=' '{print $2}' || echo "NA")

    # ribosomal RNA & mitochondrial DNA (unique & multimapped reads)
    rRNA_mtDNA_metrics="${output_dir}/rRNA_mtDNA/${sample_name}_mt_rrna_metrics.txt"
    perc_rrna=$(grep "Percentage of rRNA reads" $rRNA_mtDNA_metrics | awk -F ',' '{print $NF}' || echo "NA")
    perc_mt=$(grep "Percentage of mtDNA reads (of mapped reads)" $rRNA_mtDNA_metrics | awk -F ',' '{print $NF}' || echo "NA")
    perc_mt_mmpa=$(grep "Percentage of mtDNA in multimapped reads (primary alignment)" $rRNA_mtDNA_metrics | awk -F ',' '{print $NF}' || echo "NA")
    perc_mt_mmaa=$(grep "Percentage of mtDNA in multimapped reads (all alignments)" $rRNA_mtDNA_metrics | awk -F ',' '{print $NF}' || echo "NA")
    
    # Noise Percentage
    total_umis_all=$(grep "yesUMIs" "${solo_dir}/Features.stats" | awk '{print $2}' || echo "NA")
    total_umis_cells=$(grep "UMIs in Cells" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
    noise=$(echo "scale=4; ($total_umis_all - $total_umis_cells) / $total_umis_all" | bc || echo "NA")
    echo "done noise"

    # Percentage Intronic Reads
    exonic_sum=$(cat ${map_dir}/${sample_name}_Solo.out/Gene/raw/matrix.mtx | awk 'NR>3 {sum+=$3} END{print sum}' || echo "0")
    fullgene_sum=$(cat ${map_dir}/${sample_name}_Solo.out/GeneFull/raw/matrix.mtx | awk 'NR>3 {sum+=$3} END{print sum}' || echo "0")
    intronic_sum=$(( fullgene_sum - exonic_sum ))

    fraction_intronic=$(awk -v i="$intronic_sum" -v t="$fullgene_sum" \
                        'BEGIN{printf("%.4f", (i)/t)}' || echo "NA")
    echo "done intronic"
    
    # 10x_saturate results
    saturation_dir="${output_dir}/saturation/${sample_name}"
    saturate_07=$(cat $saturation_dir/saturation.log | grep "approximately: " | awk '{print $(NF-2) " M"}' || echo "NA")
    echo "done saturation"

    # Append results to the output file
    echo -e "STARsolo\t${sample_name}\t$n_reads\t$R1_Q30\t$R2_Q30\t\
            $n_uniquely_mapped\t$p_uniquely_mapped\t$p_multi_mapped\t$p_multi_too_many\t\
            $p_unmapped_short\t$p_unmapped_other\t\t\t$N_cells\t$UMI_cutoff\t$saturation\t\
            $saturate_07\t$noise\t$fraction_intronic\t$perc_rrna\t$perc_mt\t$perc_mt_mmaa\t$perc_mt_mmpa\t\
            $mean_reads\t$median_umis\t$median_genes\t$total_genes" >> "$output_file"
done

# ------------------------------------------------------------------
# Extract metrics from ParseBio_pipeline directory if it exists
# ------------------------------------------------------------------
if [ -d "${output_dir}/ParseBio_pipeline" ]; then
    echo "mapping dir: ${output_dir}/ParseBio_pipeline"
    
    # For each sample directory under mapping_splitpipe
    for splitpipe_dir in ${output_dir}/ParseBio_pipeline/*; do
        # Check for stats file
        report_dir="${splitpipe_dir}/all-sample/report"
        stats_file="${report_dir}/sample_all_stats.csv"
        [ -f "$stats_file" ] || continue

        # Extract fields from sample_all_stats.csv
        cDNA_Q30=$(grep "^cDNA_Q30," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        number_of_cells=$(grep "^number_of_cells," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_align_input=$(grep "^reads_align_input," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_align_multimap=$(grep "^reads_align_multimap," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_align_unique=$(grep "^reads_align_unique," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_too_many_loci=$(grep "^reads_too_many_loci," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_too_short=$(grep "^reads_too_short," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_tso_trim=$(grep "^reads_tso_trim," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        sequencing_saturation=$(grep "^sequencing_saturation," "$stats_file" | awk -F',' '{print $2}' || echo "NA")

        # Calculate percentages 
        p_multi_mapped=$(awk -v multi=$reads_align_multimap -v input=$reads_align_input 'BEGIN{if(input>0){printf("%.4f", (multi/input))}else{print("NA")}}' || echo "NA")
        p_uniquely_mapped=$(awk -v uniq=$reads_align_unique -v input=$reads_align_input 'BEGIN{if(input>0){printf("%.4f", (uniq/input))}else{print("NA")}}' || echo "NA")
        p_too_many_loci=$(awk -v val=$reads_too_many_loci -v input=$reads_align_input 'BEGIN{if(input>0){printf("%.4f", (val/input))}else{print("NA")}}' || echo "NA")
        p_too_short=$(awk -v val=$reads_too_short -v input=$reads_align_input 'BEGIN{if(input>0){printf("%.4f", (val/input))}else{print("NA")}}' || echo "NA")
        p_tso_trim=$(awk -v val=$reads_tso_trim -v input=$reads_align_input 'BEGIN{if(input>0){printf("%.4f", (val/input))}else{print("NA")}}' || echo "NA")

        # Check for agg_samp file
        agg_samp_file="${splitpipe_dir}/agg_samp_ana_summary.csv"
        [ -f "$agg_samp_file" ] || continue

        # Extract fields from agg_samp_ana_summary.csv
        mean_reads=$(grep "mean_reads_per_cell" "$agg_samp_file" | awk -F ',' '{print $NF}' || echo "NA")
        median_umis=$(grep "ref-splitpipe_median_tscp_per_cell" "$agg_samp_file" | awk -F ',' '{print $NF}' || echo "NA")
        median_genes=$(grep "ref-splitpipe_median_genes_per_cell" "$agg_samp_file" | awk -F ',' '{print $NF}' || echo "NA")

        # Append results to the output file
        echo -e "ParseBio_pipeline\t$(basename "$splitpipe_dir")\t$reads_align_input\t$cDNA_Q30\t\t$reads_align_unique\t\
                $p_uniquely_mapped\t$p_multi_mapped\t$p_too_many_loci\t$p_too_short\t$p_tso_trim\t\t\t$number_of_cells\t\
                \t$sequencing_saturation\t\
                \t\t\t\t\t\t\t\
                $mean_reads\t$median_umis\t$median_genes\t" >> "$output_file"
    done
fi

# ------------------------------------------------------------------
# Extract metrics from CellRanger_pipeline directory if it exists
# ------------------------------------------------------------------
if [ -d "${output_dir}/CellRanger_pipeline" ]; then
    echo "mapping dir: ${output_dir}/CellRanger_pipeline"
    
    # For each sample directory within the CellRanger_pipeline folder
    for cellranger_dir in ${output_dir}/CellRanger_pipeline/*; do
        # Check for metrics_summary.csv file
        metrics_file="${cellranger_dir}/outs/metrics_summary.csv"
        [ -f "$metrics_file" ] || continue

        sample_name=$(basename "$cellranger_dir" _count)
        echo "sample: $sample_name"

        # Extract fields from metrics_summary.csv
        n_cells=$(sed -n '2p' $metrics_file | sed -E 's/([0-9]),([0-9])/\1\2/g' | cut -d',' -f1 | tr -d '"' || echo "NA")
        mean_reads=$(sed -n '2p' $metrics_file | sed -E 's/([0-9]),([0-9])/\1\2/g' | cut -d',' -f2 | tr -d '"' || echo "NA")
        median_genes=$(sed -n '2p' $metrics_file | sed -E 's/([0-9]),([0-9])/\1\2/g' | cut -d',' -f3 | tr -d '"' || echo "NA")
        n_reads=$(sed -n '2p' $metrics_file | sed -E 's/([0-9]),([0-9])/\1\2/g' | cut -d',' -f4 | tr -d '"' || echo "NA")
        saturation=$(sed -n '2p' $metrics_file | sed -E 's/([0-9]),([0-9])/\1\2/g' | cut -d',' -f7 || echo "NA")
        R1_Q30=$(sed -n '2p' $metrics_file | sed -E 's/([0-9]),([0-9])/\1\2/g' | cut -d',' -f9 || echo "NA")
        p_uniquely_mapped=$(sed -n '2p' $metrics_file | sed -E 's/([0-9]),([0-9])/\1\2/g' | cut -d',' -f12 | tr -d '"' || echo "NA")
        total_genes=$(sed -n '2p' $metrics_file | sed -E 's/([0-9]),([0-9])/\1\2/g' | cut -d',' -f19 | tr -d '"' || echo "NA")
        median_umis=$(sed -n '2p' $metrics_file | sed -E 's/([0-9]),([0-9])/\1\2/g' | cut -d',' -f20 | tr -d '"' || echo "NA")

        # Append results to the output file
        echo -e "CellRanger_pipeline\t$sample_name\t$n_reads\t$R1_Q30\t\t\
                \t$p_uniquely_mapped\t\t\t\
                \t\t\t\t$n_cells\t\t$saturation\t\
                \t\t\t\t\t\t\t\
                $mean_reads\t$median_umis\t$median_genes\t$total_genes" >> "$output_file"

    done
fi

echo "Finished. Results are in $output_file"
