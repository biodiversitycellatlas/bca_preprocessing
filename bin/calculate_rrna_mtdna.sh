#!/bin/bash

# ------------------------------------------------------------------
# Define Variables
# ------------------------------------------------------------------
bam_file=$1
outfile=$2
ref_star_gtf=$3
grep_rrna=$4
mt_contig=$5


# ------------------------------------------------------------------
# Create output file and general information
# ------------------------------------------------------------------

# Create output file
echo -e "Metric,Count" > $outfile
echo -e "GTF file,${ref_star_gtf}" >> $outfile
echo -e "rRNA Grep command,${grep_rrna}" >> $outfile
echo -e "MT Contig,${mt_contig}" >> $outfile

# Total number of mapped reads (excluding unmapped reads)
mapped=$(samtools view -F 4 $bam_file | wc -l)
echo -e "Total number of mapped reads,$mapped" >> $outfile

# Total number of unmapped reads
unmapped=$(samtools view -f 4 $bam_file | wc -l)
echo -e "Total number of unmapped reads,$unmapped" >> $outfile

# Number of uniquely mapped reads
uniquely_mapped=$(samtools view -F 4 $bam_file | grep 'NH:i:1' | wc -l)
echo -e "Number of uniquely mapped reads,$uniquely_mapped" >> $outfile


# ------------------------------------------------------------------
# Calculate rRNA metrics
# ------------------------------------------------------------------
# Number of ribosomal RNA reads in uniquely mapped reads
featureCounts -t "${grep_rrna}" -a ${ref_star_gtf} -o feat_counts_rRNA.txt $bam_file
rrna=$(cat feat_counts_rRNA.txt.summary | grep -v "Assigned" | awk '{print $2}')
echo -e "Number of ribosomal RNA reads in uniquely mapped reads,$rrna" >> $outfile

# Calculate percentage of rRNA reads among uniquely mapped reads
if [ "$uniquely_mapped" -gt 0 ]; then
    perc_rrna=$(awk -v r="$rrna" -v u="$uniquely_mapped" 'BEGIN {printf "%.2f", (r/u)*100}')
else
    perc_rrna="NA"
fi
echo -e "Percentage of rRNA reads (of uniquely mapped reads),$perc_rrna" >> $outfile


# ------------------------------------------------------------------
# Calculate mtDNA metrics
# ------------------------------------------------------------------
# Number of reads mapping to mtDNA contig
mt=$(samtools view $bam_file | grep ${mt_contig} | wc -l)
echo -e "Number of reads mapping to mtDNA contig,$mt" >> $outfile

# Percentage of mtDNA reads among all mapped reads
if [ "$mapped" -gt 0 ]; then
    perc_mt=$(awk -v m="$mt" -v tot="$mapped" 'BEGIN {printf "%.2f", (m/tot)*100}')
else
    perc_mt="NA"
fi
echo -e "Percentage of mtDNA reads (of mapped reads),$perc_mt" >> $outfile

# Multimapped reads (primary alignment)
samtools view -h -F 256 $bam_file | grep -E "^\@|NH:i:[2-9]" | samtools view -b -o multimapped_primealign.bam
total_mmpa=$(samtools view multimapped_primealign.bam | wc -l)
echo -e "Total multimapped reads (primary alignment),$total_mmpa" >> $outfile
mt_multi1=$(samtools view multimapped_primealign.bam | grep ${mt_contig} | wc -l)
echo -e "mtDNA counts in Multimapped reads (primary alignment),\$mt_multi1" >> $outfile

# Percentage of mtDNA reads among multimapped reads (primary alignment)
if [ "$total_mmpa" -gt 0 ]; then
    perc_mt_mmpa=$(awk -v m="$mt_multi1" -v tot="$total_mmpa" 'BEGIN {printf "%.2f", (m/tot)*100}')
else
    perc_mt_mmpa="NA"
fi
echo -e "Percentage of mtDNA in multimapped reads (primary alignment),$perc_mt_mmpa" >> $outfile

# Multimapped reads (all alignments)
samtools view -h $bam_file | grep -E "^\@|NH:i:[2-9]" | samtools view -b -o multimapped_allalign.bam
total_mmaa=$(samtools view multimapped_allalign.bam | wc -l)
echo -e "Total multimapped reads (all alignments),$total_mmaa" >> $outfile
mt_multi2=$(samtools view multimapped_allalign.bam | grep ${mt_contig} | wc -l)
echo -e "mtDNA counts in Multimapped reads (all alignments),$mt_multi2" >> $outfile

# Percentage of mtDNA reads among multimapped reads (all alignments)
if [ "$total_mmaa" -gt 0 ]; then
    perc_mt_mmaa=$(awk -v m="$mt_multi2" -v tot="$total_mmaa" 'BEGIN {printf "%.2f", (m/tot)*100}')
else
    perc_mt_mmaa="NA"
fi
echo -e "Percentage of mtDNA in multimapped reads (all alignments),$perc_mt_mmaa" >> $outfile
