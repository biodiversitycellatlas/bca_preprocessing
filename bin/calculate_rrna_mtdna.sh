#!/bin/bash

# ------------------------------------------------------------------
# Define Variables
# ------------------------------------------------------------------
bam_file=$1
outfile=$2
ref_gtf=$3
mt_contig=$4

# ------------------------------------------------------------------
# Extract rRNA biotype labels from GTF
# ------------------------------------------------------------------
# Collects all distinct gene_biotype values containing "rRNA"
rrna_biotypes=$(grep -w 'gene_biotype' "${ref_gtf}" \
    | grep -oi 'gene_biotype "[^"]*"' \
    | grep -i 'rRNA' \
    | awk -F'"' '{print $2}' \
    | sort -u \
    | tr '\n' '|' \
    | sed 's/|$//')

# ------------------------------------------------------------------
# Create output file and general information
# ------------------------------------------------------------------
echo -e "Metric,Count" > $outfile
echo -e "GTF file,${ref_gtf}" >> $outfile
echo -e "rRNA biotypes detected,${rrna_biotypes}" >> $outfile
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
# featureCounts groups counts by gene_biotype; we sum all rRNA biotype rows.

# Uniquely mapped
featureCounts -t gene -g gene_biotype -a ${ref_gtf} -o feat_counts_rRNA.txt $bam_file
rrna=$(grep -E "^(${rrna_biotypes})" feat_counts_rRNA.txt | awk '{sum += $NF} END {print sum+0}')
echo -e "Number of ribosomal RNA reads in uniquely mapped reads,$rrna" >> $outfile

perc_rrna=$(awk -v r="$rrna" -v u="$uniquely_mapped" 'BEGIN {printf "%.4f", (r/u)}')
echo -e "Percentage of rRNA reads (of uniquely mapped reads),$perc_rrna" >> $outfile

# Multimapped primary alignments
samtools view -h -F 256 $bam_file | grep -E "^\@|NH:i:[2-9]" | samtools view -b -o multimapped_primealign.bam
total_mmpa=$(samtools view multimapped_primealign.bam | wc -l)
echo -e "Total multimapped reads (primary alignment),$total_mmpa" >> $outfile

featureCounts -t gene -g gene_biotype -a ${ref_gtf} -o feat_counts_rRNA_mmpa.txt multimapped_primealign.bam
rrna_mmpa=$(grep -E "^(${rrna_biotypes})" feat_counts_rRNA_mmpa.txt | awk '{sum += $NF} END {print sum+0}')
echo -e "rRNA counts in Multimapped reads (primary alignment),$rrna_mmpa" >> $outfile

perc_rrna_mmpa=$(awk -v r="$rrna_mmpa" -v tot="$total_mmpa" 'BEGIN {printf "%.4f", (r/tot)}')
echo -e "Percentage of rRNA in multimapped reads (primary alignment),$perc_rrna_mmpa" >> $outfile

# Multimapped all alignments
samtools view -h $bam_file | grep -E "^\@|NH:i:[2-9]" | samtools view -b -o multimapped_allalign.bam
total_mmaa=$(samtools view multimapped_allalign.bam | wc -l)
echo -e "Total multimapped reads (all alignments),$total_mmaa" >> $outfile

featureCounts -t gene -g gene_biotype -a ${ref_gtf} -o feat_counts_rRNA_mmaa.txt multimapped_allalign.bam
rrna_mmaa=$(grep -E "^(${rrna_biotypes})" feat_counts_rRNA_mmaa.txt | awk '{sum += $NF} END {print sum+0}')
echo -e "rRNA counts in Multimapped reads (all alignments),$rrna_mmaa" >> $outfile

perc_rrna_mmaa=$(awk -v r="$rrna_mmaa" -v tot="$total_mmaa" 'BEGIN {printf "%.4f", (r/tot)}')
echo -e "Percentage of rRNA in multimapped reads (all alignments),$perc_rrna_mmaa" >> $outfile

# ------------------------------------------------------------------
# Calculate mtDNA metrics
# ------------------------------------------------------------------
# Number of reads mapping to mtDNA contig
mt=$(samtools view $bam_file | grep ${mt_contig} | wc -l)
echo -e "Number of reads mapping to mtDNA contig,$mt" >> $outfile

# Percentage of mtDNA reads among all mapped reads
perc_mt=$(awk -v m="$mt" -v tot="$mapped" 'BEGIN {printf "%.4f", (m/tot)}')
echo -e "Percentage of mtDNA reads (of mapped reads),$perc_mt" >> $outfile

# Multimapped reads (primary alignment)
samtools view -h -F 256 $bam_file | grep -E "^\@|NH:i:[2-9]" | samtools view -b -o multimapped_primealign.bam
total_mmpa=$(samtools view multimapped_primealign.bam | wc -l)
echo -e "Total multimapped reads (primary alignment),$total_mmpa" >> $outfile

mt_multi1=$(samtools view multimapped_primealign.bam | grep ${mt_contig} | wc -l)
echo -e "mtDNA counts in Multimapped reads (primary alignment),$mt_multi1" >> $outfile

# Percentage of mtDNA reads among multimapped reads (primary alignment)
perc_mt_mmpa=$(awk -v m="$mt_multi1" -v tot="$total_mmpa" 'BEGIN {printf "%.4f", (m/tot)}')
echo -e "Percentage of mtDNA in multimapped reads (primary alignment),$perc_mt_mmpa" >> $outfile

# Multimapped reads (all alignments)
samtools view -h $bam_file | grep -E "^\@|NH:i:[2-9]" | samtools view -b -o multimapped_allalign.bam
total_mmaa=$(samtools view multimapped_allalign.bam | wc -l)
echo -e "Total multimapped reads (all alignments),$total_mmaa" >> $outfile

mt_multi2=$(samtools view multimapped_allalign.bam | grep ${mt_contig} | wc -l)
echo -e "mtDNA counts in Multimapped reads (all alignments),$mt_multi2" >> $outfile

# Percentage of mtDNA reads among multimapped reads (all alignments)
perc_mt_mmaa=$(awk -v m="$mt_multi2" -v tot="$total_mmaa" 'BEGIN {printf "%.4f", (m/tot)}')
echo -e "Percentage of mtDNA in multimapped reads (all alignments),$perc_mt_mmaa" >> $outfile
