#!/usr/bin/env bash
#
# Barcode check, run on every sample before mapping.
#
# The mapping-rate check (bin/mapping_rate_guard.sh) does not depend on the barcodes: with
# a whitelist, a protocol or a barcode orientation that does not match the reads, the cDNA
# maps as well as ever and the run ends hours later with an almost empty matrix. STARsolo
# writes its barcode statistics only when it finishes, so this reads the barcodes itself:
# the first --reads barcode reads are cut at the positions of the layout STAR is given,
# and each segment is matched exactly against its whitelist.
#
# The verdict is written to --out, never returned as an exit code, so a failed sample
# still hands its numbers to the pipeline (lib/BarcodeCheck.groovy stops the run or drops
# the sample):
#
#   pass     the share of reads whose barcode is on the whitelist is at least --min
#   fail     the share is below --min, or the barcode read is shorter than the layout
#   skipped  nothing to check: no whitelist, a layout that is not at fixed positions from
#            the read start, a whitelist count that differs from the barcode segments, or
#            too few reads
#
# Matching is exact on purpose. A random 16-mer matches a list of 6.8 million barcodes
# with one mismatch in about 8% of reads, but exactly in 0.2%, so a wrong and a right
# whitelist are far apart.
#
# Only bash, awk and pigz or gzip are used: the STAR container has no python.
#
# Usage:
#   barcode_check.sh --sample ID --fastq BARCODE_READ.fastq.gz --layout "STAR SOLO OPTIONS" \
#       --min PCT --out FILE [--reads N] [--min-reads N] [WHITELIST ...]
#
#   --layout   the --soloType ... options STAR gets, e.g.
#              "--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12"
#              "--soloType CB_UMI_Complex --soloCBposition 0_50_0_57 0_30_0_37 0_10_0_17 --soloUMIposition 0_0_0_9"
#   WHITELIST  one file per barcode segment, in segment order; plain or gzipped

set -uo pipefail
export LC_ALL=C

sample="" fastq="" layout="" min="" out="" reads=100000 min_reads=1000
whitelists=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)    sample="$2"; shift 2 ;;
        --fastq)     fastq="$2"; shift 2 ;;
        --layout)    layout="$2"; shift 2 ;;
        --min)       min="$2"; shift 2 ;;
        --out)       out="$2"; shift 2 ;;
        --reads)     reads="$2"; shift 2 ;;
        --min-reads) min_reads="$2"; shift 2 ;;
        -h|--help)   sed -n '2,37p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        --*)         echo "barcode_check.sh: unknown option $1" >&2; exit 2 ;;
        *)           whitelists+=("$1"); shift ;;
    esac
done
if [[ -z "$sample" || -z "$fastq" || -z "$min" || -z "$out" ]]; then
    echo "barcode_check.sh: --sample, --fastq, --min and --out are required" >&2
    exit 2
fi
for f in "$fastq" ${whitelists[@]+"${whitelists[@]}"}; do
    [[ -r "$f" ]] || { echo "barcode_check.sh: cannot read $f" >&2; exit 2; }
done

# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

# write_verdict VERDICT READS PCT_MATCH PCT_REVCOMP PCT_PER_SEGMENT READ_LENGTH LAYOUT_LENGTH NOTE
write_verdict() {
    {
        printf 'sample\tverdict\treads\tpct_match\tpct_revcomp\tpct_per_segment\tread_length\tlayout_length\tthreshold\tnote\n'
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$min" "$8"
    } >"$out"
    echo "[barcode_check] ${sample}: $1${8:+ ($8)}"
}

# skip NOTE [LAYOUT_LENGTH] -- nothing to check: the sample goes on
skip() { write_verdict skipped NA NA NA NA NA "${2:-NA}" "$1"; exit 0; }

# below A B -- floating-point A < B
below() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 < b + 0) }'; }

# Plain text of a plain or gzipped file. A gzip file starts with the bytes 1f 8b.
open_text() {
    if [[ "$(head -c 2 "$1")" == $'\x1f\x8b' ]]; then
        if command -v pigz >/dev/null 2>&1; then pigz -dc "$1"; else gzip -dc "$1"; fi
    else
        cat "$1"
    fi
}

# --------------------------------------------------------------------------
# Layout: the barcode segments as 1-based start and length, and where the layout ends
# --------------------------------------------------------------------------

solo_type="" cb_start=1 cb_len=16 umi_start=17 umi_len=10 umi_pos=""
cb_pos=()
read -r -a tok <<<"$layout"
i=0
while (( i < ${#tok[@]} )); do
    case "${tok[i]}" in
        --soloType)        solo_type="${tok[i+1]:-}"; i=$((i + 2)) ;;
        --soloCBstart)     cb_start="${tok[i+1]:-}";  i=$((i + 2)) ;;
        --soloCBlen)       cb_len="${tok[i+1]:-}";    i=$((i + 2)) ;;
        --soloUMIstart)    umi_start="${tok[i+1]:-}"; i=$((i + 2)) ;;
        --soloUMIlen)      umi_len="${tok[i+1]:-}";   i=$((i + 2)) ;;
        --soloUMIposition) umi_pos="${tok[i+1]:-}";   i=$((i + 2)) ;;
        --soloCBposition)
            i=$((i + 1))
            while (( i < ${#tok[@]} )) && [[ "${tok[i]}" != --* ]]; do cb_pos+=("${tok[i]}"); i=$((i + 1)); done ;;
        *) i=$((i + 1)) ;;
    esac
done

seg_start=() seg_len=()
layout_len=""
case "$solo_type" in
    CB_UMI_Simple)
        seg_start=("$cb_start"); seg_len=("$cb_len")
        layout_len=$(( cb_start + cb_len - 1 ))
        (( umi_start + umi_len - 1 > layout_len )) && layout_len=$(( umi_start + umi_len - 1 ))
        ;;
    CB_UMI_Complex)
        # startAnchor_startPosition_endAnchor_endPosition, positions 0-based. Anchor 0 is
        # the read start; the others (read end, adapter) move from read to read.
        layout_len=0
        for pos in ${cb_pos[@]+"${cb_pos[@]}"}; do
            IFS=_ read -r a b c d <<<"$pos"
            if [[ "$a" != 0 || "$c" != 0 || ! "$b" =~ ^[0-9]+$ || ! "$d" =~ ^[0-9]+$ ]]; then
                skip "the barcode position $pos is not counted from the read start"
            fi
            seg_start+=($(( b + 1 ))); seg_len+=($(( d - b + 1 )))
            (( d + 1 > layout_len )) && layout_len=$(( d + 1 ))
        done
        if [[ -n "$umi_pos" ]]; then
            IFS=_ read -r a b c d <<<"$umi_pos"
            if [[ "$a" == 0 && "$c" == 0 && "$d" =~ ^[0-9]+$ ]] && (( d + 1 > layout_len )); then
                layout_len=$(( d + 1 ))
            fi
        fi
        ;;
    *)
        skip "no barcode layout to check (soloType '${solo_type:-not given}')" ;;
esac
n_seg=${#seg_start[@]}
(( n_seg > 0 )) || skip "the layout names no barcode segment"

(( ${#whitelists[@]} > 0 )) || skip "no barcode whitelist" "$layout_len"
if (( ${#whitelists[@]} != n_seg )); then
    skip "${#whitelists[@]} whitelist file(s) for ${n_seg} barcode segment(s)" "$layout_len"
fi

# --------------------------------------------------------------------------
# Match
# --------------------------------------------------------------------------

# One stream for awk: the whitelists as "W <segment> <barcode>", then the sequence line
# of the first --reads records as "R <sequence>". head closes the pipe once it has its
# lines, which ends the decompression early; that exit is expected, not an error.
result=$(
    {
        s=0
        for wl in "${whitelists[@]}"; do
            s=$((s + 1))
            open_text "$wl" | tr -d '\r' | awk -v s="$s" 'NF { print "W\t" s "\t" toupper($1) }'
        done
        open_text "$fastq" 2>/dev/null | head -n $(( 4 * reads )) | awk 'NR % 4 == 2 { print "R\t" $0 }'
    } | awk -F'\t' -v starts="${seg_start[*]}" -v lens="${seg_len[*]}" '
        function revcomp(s,    i, out) {
            out = ""
            for (i = length(s); i >= 1; i--) out = out comp[substr(s, i, 1)]
            return out
        }
        BEGIN {
            n_seg = split(starts, S, " "); split(lens, L, " ")
            comp["A"] = "T"; comp["C"] = "G"; comp["G"] = "C"; comp["T"] = "A"; comp["N"] = "N"
        }
        $1 == "W" {
            if (!(($2, $3) in wl)) { wl[$2, $3] = 1; wl_n[$2]++ }
            if (!($2 in wl_len)) wl_len[$2] = length($3)
            next
        }
        $1 == "R" {
            n++
            len_count[length($2)]++
            all = 1; all_rc = 1
            for (i = 1; i <= n_seg; i++) {
                bc = substr($2, S[i], L[i])
                if ((i, bc) in wl) seg[i]++; else all = 0
                if (!((i, revcomp(bc)) in wl)) all_rc = 0
            }
            if (all) matched++
            if (all_rc) matched_rc++
        }
        END {
            # The most common read length
            for (len in len_count) if (len_count[len] > best) { best = len_count[len]; mode = len }
            per_seg = ""; sizes = ""; widths = ""
            for (i = 1; i <= n_seg; i++) {
                sep = (i > 1) ? "," : ""
                per_seg = per_seg sep sprintf("%.1f", n ? 100 * seg[i] / n : 0)
                sizes   = sizes sep (wl_n[i] + 0)
                widths  = widths sep (wl_len[i] + 0)
            }
            printf "%d\t%.1f\t%.1f\t%s\t%d\t%s\t%s\n", n,
                n ? 100 * matched / n : 0, n ? 100 * matched_rc / n : 0, per_seg, mode, sizes, widths
        }'
)
IFS=$'\t' read -r n_reads pct pct_rc per_seg read_len wl_sizes wl_widths <<<"$result"
if [[ -z "${n_reads:-}" ]]; then
    echo "barcode_check.sh: no result from matching ${fastq}" >&2
    exit 1
fi

echo "[barcode_check] ${sample}: ${n_reads} reads, read length ${read_len}, layout ${layout_len} nt in ${n_seg} barcode segment(s)"
echo "[barcode_check] ${sample}: whitelist(s) of ${wl_sizes} barcodes; ${pct}% of reads on the whitelist (per segment ${per_seg}%), ${pct_rc}% as reverse complement"

(( n_reads >= min_reads )) || skip "only ${n_reads} reads, fewer than the ${min_reads} needed" "$layout_len"

# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------

# Hints for a failed sample, most specific first
hints=()
IFS=, read -r -a seg_pct <<<"$per_seg"
IFS=, read -r -a widths  <<<"$wl_widths"
for (( s = 0; s < n_seg; s++ )); do
    if (( widths[s] != seg_len[s] )); then
        hints+=("whitelist $((s + 1)) holds ${widths[s]}-nt barcodes, but segment $((s + 1)) of the layout is ${seg_len[s]} nt")
    fi
done
if ! below "$pct_rc" "$min"; then
    hints+=("${pct_rc}% of the barcodes match as reverse complement: the whitelist is in the other orientation than the reads")
fi
if (( n_seg > 1 )); then
    low=()
    for (( s = 0; s < n_seg; s++ )); do below "${seg_pct[s]}" "$min" && low+=("$((s + 1))"); done
    if (( ${#low[@]} > 0 && ${#low[@]} < n_seg )); then
        hints+=("only segment(s) ${low[*]} fail to match (per segment: ${per_seg}%): check those whitelist files and their order")
    fi
fi

fail_banner() {
    {
        echo ""
        echo "====================== BARCODE CHECK FAILED ======================="
        echo "Sample ${sample}: $1"
        echo "Mapping and all downstream steps for this sample are cancelled."
        echo ""
        if (( ${#hints[@]} > 0 )); then
            echo "What the reads show:"
            printf '  - %s\n' "${hints[@]}"
            echo ""
        fi
        echo "Likely causes: the whitelist of another kit version, a wrong 'protocol',"
        echo "whitelist files given in the wrong order, or the cDNA and barcode reads"
        echo "swapped in the samplesheet."
        echo ""
        echo "To disable this check: skip_barcode_check = true"
        echo "To lower the threshold: min_barcode_match = <percent>"
        echo "==================================================================="
        echo ""
    } >&2
}

if (( read_len < layout_len )); then
    reason="the barcode read is ${read_len} nt, but the layout of this protocol needs ${layout_len} nt"
    hints=()
    fail_banner "$reason."
    write_verdict fail "$n_reads" "$pct" "$pct_rc" "$per_seg" "$read_len" "$layout_len" "$reason"
    exit 0
fi

if below "$pct" "$min"; then
    fail_banner "${pct}% of the first ${n_reads} barcode reads are on the whitelist, below the ${min}% threshold."
    note="${hints[0]:-}"
    write_verdict fail "$n_reads" "$pct" "$pct_rc" "$per_seg" "$read_len" "$layout_len" "$note"
    exit 0
fi

# A longer read is only a mismatch for CB_UMI_Simple, where the layout is the whole read
note=""
if [[ "$solo_type" == CB_UMI_Simple ]] && (( read_len > layout_len )); then
    note="the barcode read is ${read_len} nt, longer than the ${layout_len} nt of barcode and UMI: check the kit version"
fi
write_verdict pass "$n_reads" "$pct" "$pct_rc" "$per_seg" "$read_len" "$layout_len" "$note"
exit 0
