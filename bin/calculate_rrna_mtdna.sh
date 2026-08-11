#!/bin/bash

# ------------------------------------------------------------------
# Define Variables
# ------------------------------------------------------------------
bam_file=$1
outfile=$2
ref_gtf=$3
mt_contig=$4

# Multimappers are the reads with NH > 1. The digits have to be anchored: the
# plain pattern "NH:i:1" also matches "NH:i:10", so every read with 10 or more
# alignments would be classed as unique, and "NH:i:[2-9]" misses those same
# reads on the multimapper side.
NH_UNIQUE='NH:i:1([[:space:]]|$)'
NH_MULTI='NH:i:([2-9]|[0-9][0-9]+)([[:space:]]|$)'

# ------------------------------------------------------------------
# Detect how the GTF spells the biotype, and on which feature
# ------------------------------------------------------------------
# featureCounts needs a feature type (-t) and a grouping attribute (-g) that
# actually co-occur on the same GTF line. Ensembl uses gene_biotype on gene
# rows; GENCODE and NCBI-derived GTFs use gene_type, and some annotations carry
# the biotype only on transcript or exon rows. Probing for the combination
# keeps this working across all of them -- when the attribute is absent
# featureCounts aborts, and the rRNA numbers silently read as 0%.
feature_type=""
biotype_attr=""
for ftype in gene transcript exon; do
    for attr in gene_biotype gene_type transcript_biotype transcript_type; do
        if awk -F'\t' -v t="$ftype" -v a="$attr" '
                $0 ~ /^#/       { next }
                NF < 9          { next }
                $3 == t && index($9, a " \"") > 0 { found = 1; exit }
                END             { exit !found }
            ' "${ref_gtf}"; then
            feature_type="${ftype}"
            biotype_attr="${attr}"
            break 2
        fi
    done
done

# ------------------------------------------------------------------
# Extract rRNA biotype labels from GTF
# ------------------------------------------------------------------
# Collects all distinct biotype values containing "rRNA" (rRNA, Mt_rRNA,
# rRNA_pseudogene, ...); each becomes one row in the featureCounts output.
rrna_biotypes=""
if [ -n "${biotype_attr}" ]; then
    rrna_biotypes=$(grep -oi "${biotype_attr} \"[^\"]*\"" "${ref_gtf}" \
        | grep -i 'rRNA' \
        | awk -F'"' '{print $2}' \
        | sort -u \
        | tr '\n' '|' \
        | sed 's/|$//')
fi

if [ -z "${rrna_biotypes}" ]; then
    echo "WARNING: no rRNA biotype found in ${ref_gtf}" >&2
    echo "WARNING: biotype attribute detected: '${biotype_attr:-none}', feature type: '${feature_type:-none}'" >&2
    echo "WARNING: all rRNA metrics will be reported as N/A" >&2
fi

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
# Sum the counts of every rRNA biotype row in a featureCounts table. Rows are
# matched on the exact biotype so "rRNA" does not also pick up the separate
# "rRNA_pseudogene" row, which would double-count it. Returns N/A -- never 0 --
# when featureCounts produced no table, so a failed run is distinguishable from
# an annotation that genuinely holds no rRNA reads.
sum_rrna() {
    local table=$1
    if [ -z "${rrna_biotypes}" ] || [ ! -s "${table}" ]; then
        echo "N/A"
        return
    fi
    awk -v pat="${rrna_biotypes}" '
        BEGIN { n = split(pat, want, "|"); for (i = 1; i <= n; i++) keep[want[i]] = 1 }
        /^#/            { next }
        $1 == "Geneid"  { next }
                        { rows++ }
        ($1 in keep)    { sum += $NF }
        END             { if (rows) print sum + 0; else print "N/A" }
    ' "${table}"
}

# Percentage as a 0-1 fraction, guarding against an empty or zero denominator.
frac() {
    local num=$1 den=$2
    if [ "${num}" = "N/A" ] || [ -z "${den}" ] || [ "${den}" = "0" ]; then
        echo "N/A"
    else
        awk -v n="${num}" -v d="${den}" 'BEGIN {printf "%.4f", n / d}'
    fi
}

# Run featureCounts, tolerating failure: the caller turns a missing table into
# N/A rather than aborting the mtDNA metrics, which are computed independently.
count_biotypes() {
    local out=$1 bam=$2
    shift 2
    if [ -z "${biotype_attr}" ]; then
        return
    fi
    featureCounts -t "${feature_type}" -g "${biotype_attr}" \
        -a "${ref_gtf}" -o "${out}" "$@" "${bam}" \
        || echo "WARNING: featureCounts failed for ${bam}" >&2
}

# ------------------------------------------------------------------
# Create output file and general information
# ------------------------------------------------------------------
echo -e "Metric,Count" > $outfile
echo -e "GTF file,${ref_gtf}" >> $outfile
echo -e "Biotype attribute used,${biotype_attr:-none}" >> $outfile
echo -e "Feature type used,${feature_type:-none}" >> $outfile
echo -e "rRNA biotypes detected,${rrna_biotypes}" >> $outfile
echo -e "MT Contig,${mt_contig}" >> $outfile

# Total number of mapped reads (excluding unmapped reads)
mapped=$(samtools view -F 4 $bam_file | wc -l)
echo -e "Total number of mapped reads,$mapped" >> $outfile

# Total number of unmapped reads
unmapped=$(samtools view -f 4 $bam_file | wc -l)
echo -e "Total number of unmapped reads,$unmapped" >> $outfile

# Number of uniquely mapped reads
uniquely_mapped=$(samtools view -F 4 $bam_file | grep -Ec "${NH_UNIQUE}")
echo -e "Number of uniquely mapped reads,$uniquely_mapped" >> $outfile

# ------------------------------------------------------------------
# Split out the multimapped alignments once, for reuse below
# ------------------------------------------------------------------
samtools view -h -F 256 $bam_file | grep -E "^@|${NH_MULTI}" | samtools view -b -o multimapped_primealign.bam
total_mmpa=$(samtools view -c multimapped_primealign.bam)

samtools view -h $bam_file | grep -E "^@|${NH_MULTI}" | samtools view -b -o multimapped_allalign.bam
total_mmaa=$(samtools view -c multimapped_allalign.bam)

# ------------------------------------------------------------------
# Calculate rRNA metrics
# ------------------------------------------------------------------
# featureCounts groups counts by biotype; we sum all rRNA biotype rows.

# Uniquely mapped -- featureCounts discards NH > 1 reads by default, so this
# run over the full BAM already restricts itself to the unique alignments.
count_biotypes feat_counts_rRNA.txt "$bam_file"
rrna=$(sum_rrna feat_counts_rRNA.txt)
echo -e "Number of ribosomal RNA reads in uniquely mapped reads,$rrna" >> $outfile
echo -e "Percentage of rRNA reads (of uniquely mapped reads),$(frac "$rrna" "$uniquely_mapped")" >> $outfile

# Multimapped primary alignments -- -M is required, otherwise featureCounts
# drops every read in these BAMs (they are all multimappers by construction)
# and the counts come out as a flat 0.
echo -e "Total multimapped reads (primary alignment),$total_mmpa" >> $outfile

count_biotypes feat_counts_rRNA_mmpa.txt multimapped_primealign.bam -M --primary
rrna_mmpa=$(sum_rrna feat_counts_rRNA_mmpa.txt)
echo -e "rRNA counts in Multimapped reads (primary alignment),$rrna_mmpa" >> $outfile
echo -e "Percentage of rRNA in multimapped reads (primary alignment),$(frac "$rrna_mmpa" "$total_mmpa")" >> $outfile

# Multimapped all alignments
echo -e "Total multimapped reads (all alignments),$total_mmaa" >> $outfile

count_biotypes feat_counts_rRNA_mmaa.txt multimapped_allalign.bam -M
rrna_mmaa=$(sum_rrna feat_counts_rRNA_mmaa.txt)
echo -e "rRNA counts in Multimapped reads (all alignments),$rrna_mmaa" >> $outfile
echo -e "Percentage of rRNA in multimapped reads (all alignments),$(frac "$rrna_mmaa" "$total_mmaa")" >> $outfile

# ------------------------------------------------------------------
# Calculate mtDNA metrics
# ------------------------------------------------------------------
# Number of reads mapping to mtDNA contig
mt=$(samtools view $bam_file | grep ${mt_contig} | wc -l)
echo -e "Number of reads mapping to mtDNA contig,$mt" >> $outfile

# Percentage of mtDNA reads among all mapped reads
echo -e "Percentage of mtDNA reads (of mapped reads),$(frac "$mt" "$mapped")" >> $outfile

# Multimapped reads (primary alignment)
echo -e "Total multimapped reads (primary alignment),$total_mmpa" >> $outfile

mt_multi1=$(samtools view multimapped_primealign.bam | grep ${mt_contig} | wc -l)
echo -e "mtDNA counts in Multimapped reads (primary alignment),$mt_multi1" >> $outfile

# Percentage of mtDNA reads among multimapped reads (primary alignment)
echo -e "Percentage of mtDNA in multimapped reads (primary alignment),$(frac "$mt_multi1" "$total_mmpa")" >> $outfile

# Multimapped reads (all alignments)
echo -e "Total multimapped reads (all alignments),$total_mmaa" >> $outfile

mt_multi2=$(samtools view multimapped_allalign.bam | grep ${mt_contig} | wc -l)
echo -e "mtDNA counts in Multimapped reads (all alignments),$mt_multi2" >> $outfile

# Percentage of mtDNA reads among multimapped reads (all alignments)
echo -e "Percentage of mtDNA in multimapped reads (all alignments),$(frac "$mt_multi2" "$total_mmaa")" >> $outfile
