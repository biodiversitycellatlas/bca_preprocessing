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


#################
# # start message #
# #################
# start_epoch=`date +%s`
# echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)


####################
# define variables #
####################	

# test: 
#   resDir="/no_backup/asebe/bvanwaardenburg/data/250115_ParseBio_Nvec_Tcas_Pliv_Cele/Nvec_BCA009_BCA010"
#   workDir="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas_sep"

resDir=$1

# Change to the directory that contains the data
cd ${resDir} || { echo "Error: Could not cd to data directory"; exit 1; }

# Define output file
output_file="mapping_stats.tsv"


###############
# run command #
###############

# Print header to output file
echo -e "Directory\tSample\tN reads/sample\tN R1 >Q30\tN R2 >Q30\tN uniquely mapped reads\t% uniquely mapped reads\t% multi-mapped reads\t% multi-mapped reads: too many\t% unmapped: too short\t% unmapped: other\tExpected % Doublets\tTarget N cells\tN cells\tUMI cutoff used for cell calling\tsaturation\tReads for 0.7 saturation\tNoise (% UMIs in non-cell barcodes)\t% Intronic reads\t% rRNA in Unique reads\t%rRNA in multimappers all pos\t%rRNA in multimappers primary pos\t% mtDNA in Unique reads\t%mtDNA in multimappers all pos\t%mtDNA in multimappers primary pos\t3 most freq genes in multimappers" > "$output_file"

# For each STARsolo mapping directory
for map_dir in ${resDir}/mapping_STARsolo/*; do
    echo "mapping dir: ${map_dir}"

    # For each sample directory under the mapping directory
    for sample_dir in $map_dir/*; do
        echo "sample dir: ${sample_dir} "

        sample_name=$(basename ${sample_dir})
        config="${map_dir##*_}"
        echo "config: $config"

        LOG="$sample_dir/${sample_name}_Log.final.out"
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
        solo_gene_dir="$sample_dir/${sample_name}_Solo.out/Gene"
        solo_genefull_dir="$sample_dir/${sample_name}_Solo.out/GeneFull"

        if [ -d "$solo_genefull_dir" ]; then
            solo_dir="$solo_genefull_dir"
        elif [ -d "$solo_gene_dir" ]; then
            solo_dir="$solo_gene_dir"
        else
            solo_dir=""
        fi
        echo "solo dir: $solo_dir"

        # N cells (if barcodes.tsv is available)
        if [ -n "$solo_dir" ] && [ -f "$solo_dir/filtered/barcodes.tsv" ]; then
            echo "getting cells"
            N_cells=$(wc -l < "$solo_dir/filtered/barcodes.tsv")
            echo $N_cells
        else
            echo "error here"
            N_cells="NA"
            break
        fi

        # Extract other metrics from summary file 
        if [ -n "$solo_dir" ] && [ -f "$solo_dir/Summary.csv" ]; then
            R1_Q30=$(grep "RNA" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
            R2_Q30=$(grep "CB+UMI" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
            UMI_cutoff=$(grep "nUMImin" "$sample_dir/${sample_name}_Log.out" | awk -F ';' '{print $2}' | head -n 1 | awk -F '=' '{print $2}' || echo "NA")
            saturation=$(grep "Saturation" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
        else
            R1_Q30="NA"
            R2_Q30="NA"
            UMI_cutoff="NA"
            saturation="NA"
            echo "error summary file"
        fi

        # ------------------------------------------------------------------
        # ribosomal RNA & mitochondrial DNA (unique & multimapped reads)
        # ------------------------------------------------------------------
        # Extract mtDNA & rRNA values
        featcounts_dir="${resDir}/rRNA_mtDNA/rRNA_mtDNA/${sample_name}"
        echo "featurecounts dir: ${featcounts_dir}"

        rRNA_summary="$featcounts_dir/feat_counts_rRNA.txt.summary"
        mtDNA_summary="$featcounts_dir/feat_counts_mtDNA.txt.summary"

        if [ -f "$rRNA_summary" ] && [ -f "$mtDNA_summary" ]; then
            reads_mmpa=$(samtools view -c "$featcounts_dir/multimapped_primealign.bam")
	        reads_mmaa=$(samtools view -c "$featcounts_dir/multimapped_allalign.bam")
            rRNA_mmpa_summary="$featcounts_dir/feat_counts_rRNA_mmpa.txt.summary"
            mtDNA_mmpa_summary="$featcounts_dir/feat_counts_mtDNA_mmpa.txt.summary"
            rRNA_mmaa_summary="$featcounts_dir/feat_counts_rRNA_mmaa.txt.summary"
            mtDNA_mmaa_summary="$featcounts_dir/feat_counts_mtDNA_mmaa.txt.summary"

            # counts in uniquely mapped reads
            rRNA_assigned=$(grep "^Assigned" "$rRNA_summary" | awk '{print $2}')
            mtDNA_assigned=$(grep "^Assigned" "$mtDNA_summary" | awk '{print $2}')
            rRNA_percentage=$(awk -v sum=$n_uniquely_mapped -v assigned=$rRNA_assigned 'BEGIN{if(assigned>0){printf("%.6f", (assigned/sum))}else{print("NA")}}')
            mtDNA_percentage=$(awk -v sum=$n_uniquely_mapped -v assigned2=$mtDNA_assigned 'BEGIN{if(assigned2>0){printf("%.6f", (assigned2/sum))}else{print("NA")}}')

            # percentage of multimappers, in only the first alignment
            rRNA_mmpa_assigned=$(grep "^Assigned" "$rRNA_mmpa_summary" | awk '{print $2}')
            mtDNA_mmpa_assigned=$(grep "^Assigned" "$mtDNA_mmpa_summary" | awk '{print $2}')
            rRNA_mmpa_percentage=$(awk -v sum2=$reads_mmpa -v assigned3=$rRNA_mmpa_assigned 'BEGIN{if(assigned3>0){printf("%.6f", (assigned3/sum2))}else{print("NA")}}')
            mtDNA_mmpa_percentage=$(awk -v sum2=$reads_mmpa -v assigned4=$mtDNA_mmpa_assigned 'BEGIN{if(assigned4>0){printf("%.6f", (assigned4/sum2))}else{print("NA")}}')

            # percentage of multimappers, in all alignments (primary + secondary in featureCounts)
            rRNA_mmaa_assigned=$(grep "^Assigned" "$rRNA_mmaa_summary" | awk '{print $2}')
            mtDNA_mmaa_assigned=$(grep "^Assigned" "$mtDNA_mmaa_summary" | awk '{print $2}')
            rRNA_mmaa_percentage=$(awk -v sum3=$reads_mmaa -v assigned5=$rRNA_mmaa_assigned 'BEGIN{if(assigned5>0){printf("%.6f", (assigned5/sum3))}else{print("NA")}}')
            mtDNA_mmaa_percentage=$(awk -v sum3=$reads_mmaa -v assigned6=$mtDNA_mmaa_assigned 'BEGIN{if(assigned6>0){printf("%.6f", (assigned6/sum3))}else{print("NA")}}')
            
            # finds most frequent genes among multimappers (all positions)
            freq_mm_genes=$(head -n 3 $featcounts_dir/sorted_multimapper_gene_counts.txt | tr '\t' ':' | tr '\n' ',')

        else
            echo "Error: rRNA or mtDNA summary file missing in $dir"
            rRNA_percentage="NA"
            mtDNA_percentage="NA"
            rRNA_mmpa_percentage="NA"
            mtDNA_mmpa_percentage="NA"
            rRNA_mmaa_percentage="NA"
            mtDNA_mmaa_percentage="NA"
            freq_mm_genes="NA"
        fi
        
        # ------------------------------------------------------------------
        # Noise Percentage
        # ------------------------------------------------------------------
        total_umis_all=$(grep "yesUMIs" "${solo_dir}/Features.stats" | awk '{print $2}' || echo "NA")
        total_umis_cells=$(grep "UMIs in Cells" "$solo_dir/Summary.csv" | awk -F ',' '{print $NF}' || echo "NA")
        noise=$(echo "scale=4; ($total_umis_all - $total_umis_cells) / $total_umis_all" | bc || echo "NA")
        echo "done noise"

        # ------------------------------------------------------------------
        # Percentage Intronic Reads
        # ------------------------------------------------------------------
        exonic_sum=$(cat ${work_dir}/${sample_name}_Solo.out/Gene/raw/matrix.mtx | awk 'NR>3 {sum+=$3} END{print sum}' || echo "0")
        fullgene_sum=$(cat ${work_dir}/${sample_name}_Solo.out/GeneFull/raw/matrix.mtx | awk 'NR>3 {sum+=$3} END{print sum}' || echo "0")
        intronic_sum=$(( fullgene_sum - exonic_sum ))

        fraction_intronic=$(awk -v i="$intronic_sum" -v t="$fullgene_sum" \
                            'BEGIN{printf("%.4f", (i)/t)}' || echo "NA")
        echo "done intronic"
        
        # ------------------------------------------------------------------
        # 10x_saturate results
        # - reads needed for 0.7 saturation
        # ------------------------------------------------------------------
        saturation_dir="${resDir}/saturation/saturation/${sample_name}"
        saturate_07=$(cat $saturation_dir/saturation.log | grep "approximately: " | awk '{print $(NF-2) " M"}' || echo "NA")
        echo "done saturation"

        # ------------------------------------------------------------------
        # Append results to the output file
        # ------------------------------------------------------------------
        echo -e "$(basename "$map_dir")\t$(basename "$sample_dir")\t$n_reads\t$R1_Q30\t$R2_Q30\t$n_uniquely_mapped\t$p_uniquely_mapped\t$p_multi_mapped\t$p_multi_too_many\t$p_unmapped_short\t$p_unmapped_other\t\t\t$N_cells\t$UMI_cutoff\t$saturation\t$saturate_07\t$noise\t$fraction_intronic\t$rRNA_percentage\t$rRNA_mmaa_percentage\t$rRNA_mmpa_percentage\t$mtDNA_percentage\t$mtDNA_mmaa_percentage\t$mtDNA_mmpa_percentage\t$freq_mm_genes" >> "$output_file"
    done
done

# Handle the mapping_splitpipe directory if it exists
if [ -d "${resDir}/mapping_splitpipe" ]; then
    echo "mapping dir: ${resDir}/mapping_splitpipe"
    
    # For each sample directory under mapping_splitpipe
    for splitpipe_dir in ${resDir}/mapping_splitpipe/*; do
        # Check for stats file
        report_dir="${splitpipe_dir}/all-sample/report"
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

        echo -e "mapping_splitpipe\t$(basename "$splitpipe_dir")\t$reads_align_input\t$cDNA_Q30\t\t$reads_align_unique\t$p_uniquely_mapped\t$p_multi_mapped\t$p_too_many_loci\t$p_too_short\t$p_tso_trim\t\t\t$Nvec_number_of_cells\t\t$sequencing_saturation\t\t\t" >> "$output_file"
    done
fi

echo "Finished. Results are in $output_file"
