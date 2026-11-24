#!/usr/bin/env bash
#
# sci_cb_code.sh -- split cell barcodes as a short code for alevin-fry, and back.
#
# sci-rocket writes a 48-nt barcode read: p7, p5, ligation and RT barcode (10 nt each),
# then the 8-nt UMI. The 40-nt cell barcode is longer than salmon alevin takes (31 nt).
# ScaleBio v1.1's 37-nt read (scalebio_demux.py) holds a 29-nt one (ligation, RT, PCR),
# for which alevin-fry's generate-permit-list takes ~8 GB. So alevin-fry gets the cell's
# position in each of the barcode lists instead, as one mixed-radix number written in
# base-4 digits (A/C/G/T, most significant first), the scheme of encode_cell() in
# scalebio_demux.py. The barcodes are already corrected (by sci-rocket, or by
# scalebio_demux.py), so every read's segments are on the lists.
#
#   check   the lists can be read and their combinations fit the code
#   encode  rewrite the barcode read as code + UMI, and list the codes that occur
#   decode  codes back to the cell barcode, its segments joined with '_' as STARsolo
#           writes a CB_UMI_Complex barcode
#
# Usage:
#   sci_cb_code.sh check  --lists L1 .. Ln --segments 1-10,11-20,.. --digits 16
#   sci_cb_code.sh encode --lists L1 .. Ln --segments 1-10,11-20,.. --umi 41-48 --digits 16 \
#                         --r1 R1.fastq[.gz] --out OUT.fastq --permit codes.txt
#   sci_cb_code.sh decode --lists L1 .. Ln --digits 16 [--with-code] < codes > barcodes
#
# The lists are given in segment order, one barcode per line; repeated entries count
# once, in the order they first appear. 'encode' writes '-' as --out for stdout.
# A read whose segment is on no list, or which is too short, gets a code of N's, which
# alevin-fry does not count, as STARsolo's exact whitelist match would not; their number
# is reported. Exit codes: 0 ok, 2 usage, 3 lists that do not fit.

set -euo pipefail

usage() {
    sed -n '3,30p' "$0" | sed 's/^# \{0,1\}//'
}

die() { echo "sci_cb_code.sh: $*" >&2; exit "${2:-2}"; }

[[ $# -ge 1 ]] || { usage >&2; exit 2; }
cmd=$1; shift
case "$cmd" in
    check|encode|decode) ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown command '$cmd' (check, encode or decode)" ;;
esac

lists=(); segments=""; umi=""; digits=16; r1=""; out=""; permit=""; with_code=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --lists)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do lists+=("$1"); shift; done ;;
        --segments)  segments=$2; shift 2 ;;
        --umi)       umi=$2; shift 2 ;;
        --digits)    digits=$2; shift 2 ;;
        --r1)        r1=$2; shift 2 ;;
        --out)       out=$2; shift 2 ;;
        --permit)    permit=$2; shift 2 ;;
        --with-code) with_code=1; shift ;;
        -h|--help)   usage; exit 0 ;;
        *) die "unknown argument '$1'" ;;
    esac
done

[[ ${#lists[@]} -ge 1 ]] || die "no --lists given"
for f in "${lists[@]}"; do [[ -r "$f" ]] || die "list '$f' cannot be read"; done
[[ "$digits" =~ ^[0-9]+$ && "$digits" -ge 1 && "$digits" -le 26 ]] \
    || die "--digits must be 1-26 (4^26 is the largest power of 4 a double holds exactly)"

# The lists, one file per segment, decompressed and without carriage returns, so awk reads
# them as plain files in segment order
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
list_files=()
for i in "${!lists[@]}"; do
    zcat -f "${lists[$i]}" | tr -d '\r' > "$tmp/list_$i"
    list_files+=("$tmp/list_$i")
done

# Shared awk: read the lists (the first nseg files) into idx[s, seq] / seqs[s, i] and
# size[s], and check that their combinations fit the code
AWK_LISTS='
function load_list(s, line,   n) {
    if (line == "") return
    if (seglen[s] && length(line) != seglen[s]) {
        printf "sci_cb_code.sh: list %d has a %d-nt entry (%s); its segment is %d nt\n", s, length(line), line, seglen[s] > "/dev/stderr"
        bad_list = 1; exit 3
    }
    # size[s] + 0: a number, so the first entry is seqs[s, 0], not seqs[s, ""]
    if (!((s, line) in idx)) { n = size[s] + 0; idx[s, line] = n; seqs[s, n] = line; size[s] = n + 1 }
}
function check_lists(   s, prod, cap) {
    prod = 1
    for (s = 1; s <= nseg; s++) {
        if (!size[s]) { printf "sci_cb_code.sh: list %d is empty\n", s > "/dev/stderr"; bad_list = 1; exit 3 }
        prod *= size[s]
    }
    cap = 4 ^ digits
    if (prod > cap) {
        printf "sci_cb_code.sh: the lists make %.0f barcode combinations, more than a %d-digit code holds (%.0f)\n", prod, digits, cap > "/dev/stderr"
        bad_list = 1; exit 3
    }
    combos = prod
}
'

nseg=${#lists[@]}

if [[ "$cmd" == "check" || "$cmd" == "encode" ]]; then
    [[ -n "$segments" ]] || die "--segments is required for '$cmd'"
    IFS=',' read -r -a seg_arr <<< "$segments"
    [[ ${#seg_arr[@]} -eq $nseg ]] || die "${#seg_arr[@]} segments in --segments, but $nseg lists"
    for s in "${seg_arr[@]}"; do
        [[ "$s" =~ ^[0-9]+-[0-9]+$ ]] || die "segment '$s' is not START-END (1-based)"
    done
fi

if [[ "$cmd" == "check" ]]; then
    awk -v nseg="$nseg" -v digits="$digits" -v segments="$segments" "$AWK_LISTS"'
    BEGIN { n = split(segments, S, ","); for (s = 1; s <= n; s++) { split(S[s], ab, "-"); seglen[s] = ab[2] - ab[1] + 1 } }
    FNR == 1 { f++ }
    { load_list(f, $1) }
    END {
        if (bad_list) exit 3
        check_lists()
        if (bad_list) exit 3
        printf "sci_cb_code.sh: %d lists (", nseg
        for (s = 1; s <= nseg; s++) printf "%s%d", (s > 1 ? " x " : ""), size[s]
        printf " = %.0f combinations) fit a %d-digit code\n", combos, digits
    }' "${list_files[@]}"
    exit
fi

if [[ "$cmd" == "encode" ]]; then
    [[ -n "$umi" && "$umi" =~ ^[0-9]+-[0-9]+$ ]] || die "--umi START-END is required for 'encode'"
    [[ -n "$r1" ]] || die "--r1 is required for 'encode'"
    [[ -r "$r1" ]] || die "--r1 '$r1' cannot be read"
    [[ -n "$out" ]] || die "--out is required for 'encode'"
    [[ -n "$permit" ]] || die "--permit is required for 'encode'"
    [[ "$out" == "-" ]] && out=/dev/stdout

    # The reads come last, on stdin, after the list files
    zcat -f "$r1" | awk -v nseg="$nseg" -v digits="$digits" -v segments="$segments" -v umi="$umi" -v permit="$permit" "$AWK_LISTS"'
    function encode(v,   out_s, d, i) {
        out_s = ""
        for (i = 0; i < digits; i++) { d = v % 4; out_s = substr("ACGT", d + 1, 1) out_s; v = (v - d) / 4 }
        return out_s
    }
    function pad(s, n, c) { while (length(s) < n) s = s c; return s }
    BEGIN {
        split(segments, S, ",")
        need = 0
        for (s = 1; s <= nseg; s++) {
            split(S[s], ab, "-"); st[s] = ab[1] + 0; seglen[s] = ab[2] - ab[1] + 1
            if (ab[2] + 0 > need) need = ab[2] + 0
        }
        split(umi, U, "-"); ust = U[1] + 0; ulen = U[2] - U[1] + 1
        if (U[2] + 0 > need) need = U[2] + 0
        ncode = pad("", digits, "N"); qcode = pad("", digits, "F")
    }
    FNR == 1 { f++; if (f == nseg + 1) { check_lists(); if (bad_list) exit 3 } }
    f <= nseg { load_list(f, $1); next }
    {
        line = FNR % 4
        if (line == 1) { header = $0; next }
        if (line == 2) { seq = $0; next }
        if (line == 3) next
        qual = $0
        nreads++
        ok = (length(seq) >= need)
        v = 0
        for (s = 1; ok && s <= nseg; s++) {
            b = substr(seq, st[s], seglen[s])
            if ((s, b) in idx) v = v * size[s] + idx[s, b]; else ok = 0
        }
        if (ok) {
            if (v in code_of) code = code_of[v]
            else { code = encode(v); code_of[v] = code; ncodes++; print code > permit }
        } else { code = ncode; nbad++ }
        u = pad(substr(seq, ust, ulen), ulen, "N"); uq = pad(substr(qual, ust, ulen), ulen, "#")
        print header; print code u; print "+"; print qcode uq
    }
    END {
        if (bad_list) exit 3
        if (f <= nseg) { check_lists(); if (bad_list) exit 3 }
        close(permit)
        printf "sci_cb_code.sh: %d reads, %d cell barcodes, %d reads with a barcode on no list (code of N)\n", nreads, ncodes, nbad > "/dev/stderr"
    }' "${list_files[@]}" - > "$out"
    # An empty read file still leaves a permit list
    [[ -e "$permit" ]] || : > "$permit"
    exit
fi

# decode
awk -v nseg="$nseg" -v digits="$digits" -v with_code="$with_code" "$AWK_LISTS"'
FNR == 1 { f++; if (f == nseg + 1) { check_lists(); if (bad_list) exit 3 } }
f <= nseg { load_list(f, $1); next }
{
    code = $1
    if (code == "") next
    if (length(code) != digits || code !~ /^[ACGT]+$/) {
        printf "sci_cb_code.sh: line %d: %s is not a %d-digit code\n", FNR, code, digits > "/dev/stderr"; bad = 1; exit 1
    }
    v = 0
    for (i = 1; i <= digits; i++) v = v * 4 + index("ACGT", substr(code, i, 1)) - 1
    if (v >= combos) {
        printf "sci_cb_code.sh: line %d: %s is beyond the barcode lists\n", FNR, code > "/dev/stderr"; bad = 1; exit 1
    }
    bc = ""
    for (s = nseg; s >= 1; s--) {
        k = v % size[s]
        bc = (s < nseg ? seqs[s, k] "_" bc : seqs[s, k])
        v = (v - k) / size[s]
    }
    print (with_code ? code "\t" bc : bc)
}
END { if (bad_list) exit 3; if (bad) exit 1 }' "${list_files[@]}" -
