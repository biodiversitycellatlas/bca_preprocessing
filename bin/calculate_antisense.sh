#!/bin/bash

# ------------------------------------------------------------------
# Antisense read fraction from a STARsolo BAM
# ------------------------------------------------------------------
# STARsolo counts a read only on the strand set by --soloStrand, so reads
# landing on the opposite strand of a gene are dropped without a trace in its
# own output. featureCounts is run twice over the same BAM -- once with the
# strand STARsolo used, once with the opposite one -- and the share of
# gene-assigned reads that only the opposite strand picks up is reported.
#
# Usage: calculate_antisense.sh BAM GTF STRAND OUTFILE [THREADS]
#   STRAND is params.star_soloStrand: Forward or Reverse. Unstranded counts
#   both strands at once, so there is no antisense to separate.

bam_file=$1
ref_gtf=$2
solo_strand=$3
outfile=$4
threads=${5:-1}

# featureCounts -s: 1 = stranded, 2 = reversely stranded
case "${solo_strand}" in
    Forward) sense_s=1; antisense_s=2 ;;
    Reverse) sense_s=2; antisense_s=1 ;;
    *)
        echo "ERROR: antisense counting needs star_soloStrand Forward or Reverse, got '${solo_strand}'" >&2
        exit 2
        ;;
esac

# Exons are counted rather than gene bodies: every GTF carries exon rows,
# whereas gene rows are missing from many non-model annotations.
feature_type="exon"

# Read the Assigned count out of a featureCounts .summary file, or N/A when
# featureCounts failed and wrote none, so a failed run is not reported as 0%.
assigned() {
    local summary=$1
    if [ ! -s "${summary}" ]; then
        echo "N/A"
        return
    fi
    awk -F'\t' '$1 == "Assigned" { print $2; found = 1 } END { if (!found) print "N/A" }' "${summary}"
}

count_strand() {
    local out=$1 strand=$2
    featureCounts -T "${threads}" -a "${ref_gtf}" -t "${feature_type}" -g gene_id \
        -s "${strand}" -o "${out}" "${bam_file}" \
        || echo "WARNING: featureCounts -s ${strand} failed for ${bam_file}" >&2
}

count_strand feat_counts_sense.txt     "${sense_s}"
count_strand feat_counts_antisense.txt "${antisense_s}"

sense=$(assigned feat_counts_sense.txt.summary)
antisense=$(assigned feat_counts_antisense.txt.summary)

# Fraction as 0-1, against the reads assigned on either strand. A read
# overlapping genes on both strands is assigned in both runs, so it sits in
# both numerator and denominator.
frac="N/A"
if [ "${sense}" != "N/A" ] && [ "${antisense}" != "N/A" ]; then
    frac=$(awk -v s="${sense}" -v a="${antisense}" \
        'BEGIN { t = s + a; if (t > 0) printf "%.4f", a / t; else print "N/A" }')
fi

echo -e "Metric,Count" > "${outfile}"
echo -e "GTF file,${ref_gtf}" >> "${outfile}"
echo -e "STARsolo strand,${solo_strand}" >> "${outfile}"
echo -e "Feature type used,${feature_type}" >> "${outfile}"
echo -e "Reads assigned sense,${sense}" >> "${outfile}"
echo -e "Reads assigned antisense,${antisense}" >> "${outfile}"
echo -e "Percentage of antisense reads (of reads assigned to genes),${frac}" >> "${outfile}"
