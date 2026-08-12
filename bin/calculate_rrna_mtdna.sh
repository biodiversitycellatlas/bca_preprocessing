#!/bin/bash

# ------------------------------------------------------------------
# Define Variables
# ------------------------------------------------------------------
# Optional leading flags, before the positional arguments.
rrna_gtf=""
while [ $# -gt 0 ]; do
    case "$1" in
        --rrna-gtf) rrna_gtf="$2"; shift 2 ;;
        --)         shift; break ;;
        -*)         echo "ERROR: unknown option $1" >&2; exit 2 ;;
        *)          break ;;
    esac
done

bam_file=$1
outfile=$2
ref_gtf=$3
shift 3
# Every remaining argument is a mitochondrial contig name. They arrive either as
# one quoted argument ("chrM M MT") or as separate words, depending on how the
# caller interpolates params.mt_contig, so both forms are flattened into one
# list here. Reading only $4 silently dropped all but the first name.
read -r -a mt_contigs <<< "$*"

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
#
# The pair annotating the MOST rows wins, rather than the first one occurring
# anywhere. References here are routinely concatenations of ref_gtf and
# ref_gtf_addfeature, and the two halves need not spell the biotype the same
# way: a handful of hand-written spike-in lines using gene_biotype must not
# outvote an entire GENCODE annotation using gene_type, or featureCounts skips
# every gene of the base annotation and the rRNA counts collapse to N/A. Ties
# keep the original order of preference.
read -r feature_type biotype_attr biotype_rows feature_rows <<< "$(awk -F'\t' '
    BEGIN {
        nf = split("gene transcript exon", ftypes, " ")
        na = split("gene_biotype gene_type transcript_biotype transcript_type", attrs, " ")
    }
    /^#/   { next }
    NF < 9 { next }
    {
        for (i = 1; i <= nf; i++) {
            if ($3 != ftypes[i]) continue
            rows[ftypes[i]]++
            for (j = 1; j <= na; j++)
                if (index($9, attrs[j] " \"") > 0) count[ftypes[i] SUBSEP attrs[j]]++
        }
    }
    END {
        best_f = ""; best_a = ""; max = 0
        for (i = 1; i <= nf; i++)
            for (j = 1; j <= na; j++) {
                c = count[ftypes[i] SUBSEP attrs[j]]
                if (c > max) { max = c; best_f = ftypes[i]; best_a = attrs[j] }
            }
        if (max > 0) print best_f, best_a, max, rows[best_f]
    }
' "${ref_gtf}")"

# Partial coverage means part of the annotation uses a different convention and
# is invisible to featureCounts under the attribute chosen above. No single
# attribute can cover such a reference, so this is reported rather than fixed.
if [ -n "${biotype_attr}" ] && [ "${biotype_rows:-0}" -lt "${feature_rows:-0}" ]; then
    echo "WARNING: ${biotype_attr} annotates only ${biotype_rows} of ${feature_rows} ${feature_type} rows in ${ref_gtf}" >&2
    echo "WARNING: rows spelling the biotype differently are not counted -- typical when" >&2
    echo "WARNING: ref_gtf and ref_gtf_addfeature follow different GTF conventions" >&2
fi

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
fi

# ------------------------------------------------------------------
# Reference sequences present in the BAM
# ------------------------------------------------------------------
# Needed twice below: to size the added rRNA contigs, and to check that the
# mitochondrial names actually exist in the reference used for mapping.
sq_info=$(samtools view -H "${bam_file}" \
    | awk -F'\t' '/^@SQ/ {
        name = ""; len = ""
        for (i = 2; i <= NF; i++) {
            if ($i ~ /^SN:/) name = substr($i, 4)
            if ($i ~ /^LN:/) len  = substr($i, 4)
        }
        if (name != "") print name "\t" len
    }')

# ------------------------------------------------------------------
# Assemble the rRNA regions to count over
# ------------------------------------------------------------------
# featureCounts is handed one meta-feature, "rRNA", built from two sources:
#
#   1. every feature of the main annotation whose biotype names rRNA, under any
#      of the known spellings rather than only the one chosen above -- a
#      concatenated reference may legitimately use both;
#   2. the whole of every contig introduced by ref_gtf_addfeature, when one was
#      given. Those files are the added rRNA reference (ref_fasta_addfeature
#      carries the sequences), so every read aligning there is an rRNA read,
#      whether or not the added GTF spells out a biotype at all -- it commonly
#      does not, which is why counting by biotype alone reported none of them.
#
# Regions overlapping inside a single meta-feature are merged by featureCounts,
# so a read is still counted once when both sources cover it.
rrna_saf="rrna_regions.saf"
printf 'GeneID\tChr\tStart\tEnd\tStrand\n' > "${rrna_saf}"

if [ -n "${feature_type}" ]; then
    awk -F'\t' -v ftype="${feature_type}" '
        /^#/        { next }
        NF < 9      { next }
        $3 != ftype { next }
        {
            s = $9
            while (match(s, /(gene_biotype|gene_type|transcript_biotype|transcript_type) "[^"]*"/)) {
                if (tolower(substr(s, RSTART, RLENGTH)) ~ /rrna/) {
                    print "rRNA\t" $1 "\t" $4 "\t" $5 "\t" ($7 == "-" ? "-" : "+")
                    break
                }
                s = substr(s, RSTART + RLENGTH)
            }
        }
    ' "${ref_gtf}" >> "${rrna_saf}"
fi

added_present=""
added_missing=""
if [ -n "${rrna_gtf}" ]; then
    if [ ! -s "${rrna_gtf}" ]; then
        echo "WARNING: added rRNA annotation is empty or unreadable: ${rrna_gtf}" >&2
    else
        for contig in $(awk -F'\t' '!/^#/ && NF >= 9 { print $1 }' "${rrna_gtf}" | sort -u); do
            len=$(printf '%s\n' "${sq_info}" | awk -F'\t' -v c="${contig}" '$1 == c { print $2; exit }')
            if [ -n "${len}" ]; then
                printf 'rRNA\t%s\t1\t%s\t+\n' "${contig}" "${len}" >> "${rrna_saf}"
                added_present="${added_present} ${contig}"
            else
                added_missing="${added_missing} ${contig}"
            fi
        done
        if [ -n "${added_missing}" ]; then
            echo "WARNING: contig(s) of the added rRNA reference are absent from the BAM header:${added_missing}" >&2
            echo "WARNING: the reads were mapped against a reference built without ref_fasta_addfeature" >&2
        fi
    fi
fi

rrna_regions=$(( $(wc -l < "${rrna_saf}") - 1 ))
if [ "${rrna_regions}" -eq 0 ]; then
    echo "WARNING: no rRNA regions to count: no rRNA biotype in ${ref_gtf} and no usable" >&2
    echo "WARNING: ref_gtf_addfeature contigs; all rRNA metrics will be reported as N/A" >&2
fi

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
# Read the rRNA row out of a featureCounts table. All rRNA regions share one
# meta-feature, so there is a single row to read. Returns N/A -- never 0 -- when
# there was nothing to count or featureCounts produced no table, so a failed run
# stays distinguishable from a library that genuinely holds no rRNA reads.
sum_rrna() {
    local table=$1
    if [ "${rrna_regions}" -eq 0 ] || [ ! -s "${table}" ]; then
        echo "N/A"
        return
    fi
    awk '
        /^#/           { next }
        $1 == "Geneid" { next }
        $1 == "rRNA"   { sum += $NF; found = 1 }
        END            { if (found) print sum + 0; else print "N/A" }
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

# Run featureCounts over the assembled rRNA regions, tolerating failure: the
# caller turns a missing table into N/A rather than aborting the mtDNA metrics,
# which are computed independently.
count_rrna() {
    local out=$1 bam=$2
    shift 2
    if [ "${rrna_regions}" -eq 0 ]; then
        return
    fi
    featureCounts -F SAF -a "${rrna_saf}" -o "${out}" "$@" "${bam}" \
        || echo "WARNING: featureCounts failed for ${bam}" >&2
}

# Count the mapped alignments whose reference name is one of the mitochondrial
# contigs. The comparison is on the RNAME column and on the whole name: a grep
# over the entire SAM line also matches the contig name inside read names
# (Illumina flowcell IDs such as HMTVFDSX3 contain "MT"), inside CIGAR strings
# ("90M" contains "M"), and inside longer contig names such as chrM_alt, which
# inflates the count -- to 100% of reads in the worst case.
count_mt() {
    local bam=$1
    samtools view -F 4 "${bam}" \
        | awk -F'\t' -v names="${mt_contigs[*]}" '
            BEGIN { n = split(names, want, " "); for (i = 1; i <= n; i++) keep[want[i]] = 1 }
            ($3 in keep)  { hits++ }
            END           { print hits + 0 }
        '
}

# ------------------------------------------------------------------
# Create output file and general information
# ------------------------------------------------------------------
echo -e "Metric,Count" > $outfile
echo -e "GTF file,${ref_gtf}" >> $outfile
echo -e "Biotype attribute used,${biotype_attr:-none}" >> $outfile
echo -e "Feature type used,${feature_type:-none}" >> $outfile
echo -e "rRNA biotypes detected,${rrna_biotypes}" >> $outfile
echo -e "rRNA reference contigs (added),${added_present# }" >> $outfile
echo -e "rRNA regions counted,${rrna_regions}" >> $outfile
echo -e "MT Contig,${mt_contigs[*]}" >> $outfile

# A contig name that is not in the BAM header can only ever yield 0 reads, which
# is indistinguishable from a genuinely mitochondria-free library. Name the ones
# that do not exist instead of reporting a plausible-looking 0%.
if [ "${#mt_contigs[@]}" -eq 0 ]; then
    echo "WARNING: no mitochondrial contig given; all mtDNA metrics will be 0" >&2
else
    missing=""
    for contig in "${mt_contigs[@]}"; do
        printf '%s\n' "${sq_info}" | cut -f1 | grep -qxF "${contig}" || missing="${missing} ${contig}"
    done
    if [ -n "${missing}" ]; then
        echo "WARNING: mitochondrial contig(s) absent from the BAM header:${missing}" >&2
        echo "WARNING: check mt_contig against the reference used for mapping" >&2
    fi
fi

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
count_rrna feat_counts_rRNA.txt "$bam_file"
rrna=$(sum_rrna feat_counts_rRNA.txt)
echo -e "Number of ribosomal RNA reads in uniquely mapped reads,$rrna" >> $outfile
echo -e "Percentage of rRNA reads (of uniquely mapped reads),$(frac "$rrna" "$uniquely_mapped")" >> $outfile

# Multimapped primary alignments -- -M is required, otherwise featureCounts
# drops every read in these BAMs (they are all multimappers by construction)
# and the counts come out as a flat 0.
echo -e "Total multimapped reads (primary alignment),$total_mmpa" >> $outfile

count_rrna feat_counts_rRNA_mmpa.txt multimapped_primealign.bam -M --primary
rrna_mmpa=$(sum_rrna feat_counts_rRNA_mmpa.txt)
echo -e "rRNA counts in Multimapped reads (primary alignment),$rrna_mmpa" >> $outfile
echo -e "Percentage of rRNA in multimapped reads (primary alignment),$(frac "$rrna_mmpa" "$total_mmpa")" >> $outfile

# Multimapped all alignments
echo -e "Total multimapped reads (all alignments),$total_mmaa" >> $outfile

count_rrna feat_counts_rRNA_mmaa.txt multimapped_allalign.bam -M
rrna_mmaa=$(sum_rrna feat_counts_rRNA_mmaa.txt)
echo -e "rRNA counts in Multimapped reads (all alignments),$rrna_mmaa" >> $outfile
echo -e "Percentage of rRNA in multimapped reads (all alignments),$(frac "$rrna_mmaa" "$total_mmaa")" >> $outfile

# ------------------------------------------------------------------
# Calculate mtDNA metrics
# ------------------------------------------------------------------
# Number of reads mapping to mtDNA contig. Unmapped records are excluded here as
# they are from the denominator: --outSAMunmapped Within keeps them in the BAM.
mt=$(count_mt "$bam_file")
echo -e "Number of reads mapping to mtDNA contig,$mt" >> $outfile

# Percentage of mtDNA reads among all mapped reads
echo -e "Percentage of mtDNA reads (of mapped reads),$(frac "$mt" "$mapped")" >> $outfile

# Multimapped reads (primary alignment)
echo -e "Total multimapped reads (primary alignment),$total_mmpa" >> $outfile

mt_multi1=$(count_mt multimapped_primealign.bam)
echo -e "mtDNA counts in Multimapped reads (primary alignment),$mt_multi1" >> $outfile

# Percentage of mtDNA reads among multimapped reads (primary alignment)
echo -e "Percentage of mtDNA in multimapped reads (primary alignment),$(frac "$mt_multi1" "$total_mmpa")" >> $outfile

# Multimapped reads (all alignments)
echo -e "Total multimapped reads (all alignments),$total_mmaa" >> $outfile

mt_multi2=$(count_mt multimapped_allalign.bam)
echo -e "mtDNA counts in Multimapped reads (all alignments),$mt_multi2" >> $outfile

# Percentage of mtDNA reads among multimapped reads (all alignments)
echo -e "Percentage of mtDNA in multimapped reads (all alignments),$(frac "$mt_multi2" "$total_mmaa")" >> $outfile
