#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --qos=vshort
#SBATCH --mem=10G
#SBATCH --job-name mapping_stats

echo -e "\n\n===============  MAPPING STATISTICS  ==============="
resDir="/users/asebe/bvanwaardenburg/git/data/241204_ParseBio_Nvec_Tcas_nuclei/Nvec"

# Change to the directory that contains the data
cd ${resDir} || { echo "Error: Could not cd to script directory"; exit 1; }

# Define output file
output_file="${resDir}/star_mapping_stats.tsv"

# Print header to output file
echo -e "Directory\tSample\tN reads/sample\tN R1 >Q30\tN R2 >Q30\tN uniquely mapped reads\t% uniquely mapped reads\t% multi-mapped reads\t% multi-mapped reads: too many\t% unmapped: too short\t% unmapped: other\tExpected % Doublets\tTarget N cells\tN cells\tUMI cutoff used for cell calling\tsaturation\tNoise (% UMIs in non-cell barcodes)\tMedian % rRNA\tMedian % mtDNA" > "$output_file"

# For each STARsolo mapping directory
for map_dir in ${resDir}/mapping_STARsolo_*; do
    echo "mapping dir: ${map_dir}"

    # For each sample directory under the mapping directory
    for sample_dir in $map_dir/*; do
        echo "sample dir: ${sample_dir} "

        LOG="$sample_dir/Log.final.out"
        [ -f "$LOG" ] || continue

        # Extract from STAR Log
        n_reads=$(grep "Number of input reads" "$LOG" | awk '{print $NF}')
        n_uniquely_mapped=$(grep "Uniquely mapped reads number" "$LOG" | awk '{print $NF}')
        p_uniquely_mapped=$(grep "Uniquely mapped reads %" "$LOG" | awk '{print $NF}')
        p_multi_mapped=$(grep "% of reads mapped to multiple loci" "$LOG" | awk '{print $NF}')
        p_multi_too_many=$(grep "% of reads mapped to too many loci" "$LOG" | awk '{print $NF}')
        p_unmapped_short=$(grep "% of reads unmapped: too short" "$LOG" | awk '{print $NF}')
        p_unmapped_other=$(grep "% of reads unmapped: other" "$LOG" | awk '{print $NF}')

        # Check Gene and GeneFull directories in Solo.out
        solo_gene_dir="$sample_dir/Solo.out/Gene"
        solo_genefull_dir="$sample_dir/Solo.out/GeneFull"

        if [ -d "$solo_genefull_dir" ]; then
            solo_dir="$solo_genefull_dir"
        elif [ -d "$solo_gene_dir" ]; then
            solo_dir="$solo_gene_dir"
        else
            solo_dir=""
        fi

        # N cells (if barcodes.tsv is available)
        if [ -n "$solo_dir" ] && [ -f "$solo_dir/filtered/barcodes.tsv" ]; then
            N_cells=$(wc -l < "$solo_dir/filtered/barcodes.tsv")
        else
            N_cells="NA"
        fi

        # Extract other metrics from summary file 
        if [ -n "$solo_dir" ] && [ -f "$solo_dir/Summary.csv" ]; then
            R1_Q30=$(grep "RNA" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
            R2_Q30=$(grep "CB+UMI" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
            UMI_cutoff=$(grep "nUMImin" "$sample_dir/Log.out" | awk -F ';' '{print $2}' | head -n 1 | awk -F '=' '{print $2}' || echo "NA")
            saturation=$(grep "Saturation" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
        else
            R1_Q30="NA"
            R2_Q30="NA"
            UMI_cutoff="NA"
            saturation="NA"
        fi

        # Extract mtDNA & rRNA values
        mapName=$(basename "$sample_dir")
        pick_config=$(echo "$map_dir" | grep -oE '_(CR|N)$')

        rRNA_summary="$resDir/rRNA_mtDNA/${mapName}${pick_config}/feat_counts_rRNA.txt.summary"
        mtDNA_summary="$resDir/rRNA_mtDNA/${mapName}${pick_config}/feat_counts_mtDNA.txt.summary"

        if [ -f "$rRNA_summary" ] && [ -f "$mtDNA_summary" ]; then
            rRNA_assigned=$(grep "^Assigned" "$rRNA_summary" | awk '{print $2}')
            mtDNA_assigned=$(grep "^Assigned" "$mtDNA_summary" | awk '{print $2}')

            rRNA_sum=$(awk '{sum+=$2} END {print sum}' "$rRNA_summary")
            mtDNA_sum=$(awk '{sum2+=$2} END {print sum2}' "$mtDNA_summary")

            rRNA_percentage=$(awk -v sum=$rRNA_sum -v assigned=$rRNA_assigned 'BEGIN{if(assigned>0){printf("%.4f", (assigned/sum))}else{print("NA")}}')
            mtDNA_percentage=$(awk -v sum2=$mtDNA_sum -v assigned2=$mtDNA_assigned 'BEGIN{if(assigned2>0){printf("%.4f", (assigned2/sum2))}else{print("NA")}}')

        else
            echo "Error: rRNA or mtDNA summary file missing in $dir"
            rRNA_percentage="NA"
            mtDNA_percentage="NA"
        fi

        # Append results to the output file
        echo -e "$(basename "$map_dir")\t$(basename "$sample_dir")\t$n_reads\t$R1_Q30\t$R2_Q30\t$n_uniquely_mapped\t$p_uniquely_mapped\t$p_multi_mapped\t$p_multi_too_many\t$p_unmapped_short\t$p_unmapped_other\t\t\t$N_cells\t$UMI_cutoff\t$saturation\t\t$rRNA_percentage\t$mtDNA_percentage" >> "$output_file"
    done
done

# Handle the mapping_parsepipe directory if it exists
if [ -d "${resDir}/mapping_parsepipe" ]; then
    echo "mapping dir: ${resDir}/mapping_parsepipe"
    
    # For each sample directory under mapping_parsepipe
    for parsepipe_dir in ${resDir}/mapping_parsepipe/*; do
        # Check for stats file
        report_dir="${parsepipe_dir}/all-sample/report"
        stats_file="${report_dir}/sample_all_stats.csv"
        [ -f "$stats_file" ] || continue

        # Extract fields from sample_all_stats.csv
        cDNA_Q30=$(grep "^cDNA_Q30," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        Nvec_number_of_cells=$(grep "^number_of_cells," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_align_input=$(grep "^reads_align_input," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_align_multimap=$(grep "^reads_align_multimap," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_align_unique=$(grep "^reads_align_unique," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_too_many_loci=$(grep "^reads_too_many_loci," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_too_short=$(grep "^reads_too_short," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        reads_tso_trim=$(grep "^reads_tso_trim," "$stats_file" | awk -F',' '{print $2}' || echo "NA")
        sequencing_saturation=$(grep "^sequencing_saturation," "$stats_file" | awk -F',' '{print $2}' || echo "NA")

        # Calculate percentages if reads_align_input is numeric and > 0
        if [[ $reads_align_input =~ ^[0-9]+$ && $reads_align_input -gt 0 ]]; then
            p_multi_mapped=$(awk -v multi=$reads_align_multimap -v input=$reads_align_input 'BEGIN{if(input>0){printf("%.4f", (multi/input))}else{print("NA")}}')
            p_uniquely_mapped=$(awk -v uniq=$reads_align_unique -v input=$reads_align_input 'BEGIN{if(input>0){printf("%.4f", (uniq/input))}else{print("NA")}}')
            p_too_many_loci=$(awk -v val=$reads_too_many_loci -v input=$reads_align_input 'BEGIN{if(input>0){printf("%.4f", (val/input))}else{print("NA")}}')
            p_too_short=$(awk -v val=$reads_too_short -v input=$reads_align_input 'BEGIN{if(input>0){printf("%.4f", (val/input))}else{print("NA")}}')
            p_tso_trim=$(awk -v val=$reads_tso_trim -v input=$reads_align_input 'BEGIN{if(input>0){printf("%.4f", (val/input))}else{print("NA")}}')
        else
            p_multi_mapped="NA"
            p_uniquely_mapped="NA"
            p_too_many_loci="NA"
            p_too_short="NA"
            p_tso_trim="NA"
        fi

        echo -e "mapping_parsepipe\t$(basename "$parsepipe_dir")\t$reads_align_input\t$cDNA_Q30\t\t$reads_align_unique\t$p_uniquely_mapped\t$p_multi_mapped\t$p_too_many_loci\t$p_too_short\t$p_tso_trim\t\t\t$Nvec_number_of_cells\t\t$sequencing_saturation\t\t" >> "$output_file"
    done
fi

echo "Finished. Results are in $output_file"
