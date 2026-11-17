#!/usr/bin/env bash
# description: Check the barcode check that stops samples whose barcodes miss the whitelist
#
# BARCODE_CHECK runs bin/barcode_check.sh on every sample before mapping. It cuts the
# first barcode reads at the positions of the STARsolo layout, matches each segment
# exactly against its whitelist, and writes a verdict (pass, fail or skipped) to a TSV.
#
# The reads and whitelists are generated here, so no sequencing data, mapper or cluster
# is needed.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

CHECK="$PROJECT_ROOT/bin/barcode_check.sh"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh barcode_check [-- OPTIONS]
       tests/checks/barcode_check.sh [OPTIONS]

Validate bin/barcode_check.sh against generated barcode reads and whitelists.

Options:
  --keep          Keep the generated fixtures for inspection.
  -h, --help      Show this message.

Cases:
  simple_pass        one segment, every barcode on the whitelist
  simple_fail        one segment, another whitelist: fail, with the banner, exit 0
  threshold          a 15% share passes at --min 10 and fails at --min 20
  revcomp_hint       a reverse-complemented whitelist fails and is named as such
  complex_pass       three segments, each on its own whitelist
  complex_segment    one wrong whitelist of three: the per-segment share names it
  whitelist_width    a whitelist of another barcode length is named
  no_whitelist       no whitelist: skipped
  anchored_layout    positions counted from an adapter: skipped
  whitelist_count    one whitelist for three segments: skipped
  few_reads          fewer reads than --min-reads: skipped
  short_read         a barcode read shorter than the layout fails
  long_read          a barcode read longer than barcode + UMI passes with a note
  plain_input        an uncompressed FASTQ and a gzipped whitelist are read
  first_reads_only   only the first --reads records are matched
  missing_file       a file that cannot be read is a usage error (exit 2)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep)    KEEP=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

log_header "barcode_check"

WORKDIR="$BCA_TEST_LOGDIR/barcode_check"
rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

# make_whitelist FILE N LENGTH SEED -- N random barcodes, one per line
make_whitelist() {
    awk -v n="$2" -v len="$3" -v seed="$4" 'BEGIN {
        srand(seed); split("A C G T", B, " ")
        for (i = 0; i < n; i++) {
            s = ""; for (k = 0; k < len; k++) s = s B[int(rand() * 4) + 1]
            print s
        }
    }' >"$1"
}

# revcomp_list IN OUT -- every barcode reverse-complemented
revcomp_list() {
    awk 'BEGIN { c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" }
         { s = ""; for (i = length($1); i >= 1; i--) s = s c[substr($1, i, 1)]; print s }' "$1" >"$2"
}

# make_reads N SEED SPEC WHITELIST... -- FASTQ records on stdout.
# SPEC lists the pieces of the read: Wk is a barcode drawn from the k-th whitelist,
# Un is n random bases (a UMI or a linker).
make_reads() {
    local n="$1" seed="$2" spec="$3"; shift 3
    awk -v n="$n" -v seed="$seed" -v spec="$spec" '
        FNR == 1 { f++ }
        { wl[f, FNR] = $1; wn[f] = FNR }
        END {
            srand(seed); split("A C G T", B, " "); np = split(spec, P, ",")
            for (r = 1; r <= n; r++) {
                s = ""
                for (p = 1; p <= np; p++) {
                    v = substr(P[p], 2) + 0
                    if (substr(P[p], 1, 1) == "W") s = s wl[v, int(rand() * wn[v]) + 1]
                    else for (k = 0; k < v; k++) s = s B[int(rand() * 4) + 1]
                }
                q = s; gsub(/./, "I", q)
                printf "@read%d_%d\n%s\n+\n%s\n", seed, r, s, q
            }
        }' "$@"
}

SIMPLE="--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12"
# bc1 (8) + linker (4) + bc2 (8) + linker (4) + bc3 (8) + UMI (8)
COMPLEX="--soloType CB_UMI_Complex --soloCBposition 0_0_0_7 0_12_0_19 0_24_0_31 --soloUMIposition 0_32_0_39"

make_whitelist "$WORKDIR/wl16.txt"   500 16 1
make_whitelist "$WORKDIR/other16.txt" 500 16 2
make_whitelist "$WORKDIR/wl14.txt"   500 14 3
make_whitelist "$WORKDIR/seg1.txt"    96  8 4
make_whitelist "$WORKDIR/seg2.txt"    96  8 5
make_whitelist "$WORKDIR/seg3.txt"    96  8 6
make_whitelist "$WORKDIR/segX.txt"    96  8 7
revcomp_list "$WORKDIR/wl16.txt" "$WORKDIR/wl16_rc.txt"
gzip -c "$WORKDIR/wl16.txt" >"$WORKDIR/wl16.txt.gz"

make_reads 2000 11 W1,U12 "$WORKDIR/wl16.txt" | gzip >"$WORKDIR/simple.fastq.gz"
make_reads 2000 11 W1,U12 "$WORKDIR/wl16.txt" >"$WORKDIR/simple.fastq"
{ make_reads 300 12 W1,U12 "$WORKDIR/wl16.txt"; make_reads 1700 13 W1,U12 "$WORKDIR/other16.txt"; } \
    | gzip >"$WORKDIR/mixed.fastq.gz"
{ make_reads 1500 14 W1,U12 "$WORKDIR/wl16.txt"; make_reads 5000 15 W1,U12 "$WORKDIR/other16.txt"; } \
    | gzip >"$WORKDIR/head.fastq.gz"
make_reads 2000 16 W1,U4,W2,U4,W3,U8 "$WORKDIR/seg1.txt" "$WORKDIR/seg2.txt" "$WORKDIR/seg3.txt" \
    | gzip >"$WORKDIR/complex.fastq.gz"
make_reads 200  17 W1,U12 "$WORKDIR/wl16.txt" | gzip >"$WORKDIR/few.fastq.gz"
make_reads 2000 18 W1,U10 "$WORKDIR/wl16.txt" | gzip >"$WORKDIR/short.fastq.gz"
make_reads 2000 19 W1,U14 "$WORKDIR/wl16.txt" | gzip >"$WORKDIR/long.fastq.gz"

# run_check NAME FASTQ LAYOUT MIN [options and whitelists...]
# Sets RC, and TSV / LOG to the verdict file and the captured output of the case.
run_check() {
    local name="$1" fastq="$2" layout="$3" min="$4"; shift 4
    TSV="$WORKDIR/$name.tsv"; LOG="$WORKDIR/$name.log"
    RC=0
    bash "$CHECK" --sample S1 --fastq "$fastq" --layout "$layout" --min "$min" --out "$TSV" "$@" >"$LOG" 2>&1 || RC=$?
}

# field COLUMN -- that column of the verdict row
field() {
    awk -F'\t' -v col="$1" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == col) c = i } NR == 2 { print $c }' "$TSV"
}

expect() {  # expect NAME CONDITION_RESULT DETAIL
    if [[ "$2" -eq 0 ]]; then record PASS "$1"; else record FAIL "$1" "$3"; fi
}

summary() { echo "rc=$RC $(sed -n 2p "$TSV" 2>/dev/null | tr '\t' '|')"; }

# --------------------------------------------------------------------------
# One barcode segment
# --------------------------------------------------------------------------

log_step "one segment"

run_check simple_pass "$WORKDIR/simple.fastq.gz" "$SIMPLE" 10 "$WORKDIR/wl16.txt"
ok=0
[[ $RC -eq 0 && "$(field verdict)" == pass && "$(field pct_match)" == 100.0 ]] || ok=1
[[ "$(field reads)" == 2000 && "$(field read_length)" == 28 && "$(field layout_length)" == 28 ]] || ok=1
expect simple_pass $ok "$(summary)"

run_check simple_fail "$WORKDIR/simple.fastq.gz" "$SIMPLE" 10 "$WORKDIR/other16.txt"
ok=0
[[ $RC -eq 0 && "$(field verdict)" == fail && "$(field pct_match)" == 0.0 ]] || ok=1
grep -q "BARCODE CHECK FAILED" "$LOG" || ok=1
grep -q "skip_barcode_check" "$LOG" || ok=1
expect simple_fail $ok "$(summary)"

run_check threshold_10 "$WORKDIR/mixed.fastq.gz" "$SIMPLE" 10 "$WORKDIR/wl16.txt"
ok=0
[[ "$(field verdict)" == pass && "$(field pct_match)" == 15.0 ]] || ok=1
run_check threshold_20 "$WORKDIR/mixed.fastq.gz" "$SIMPLE" 20 "$WORKDIR/wl16.txt"
[[ "$(field verdict)" == fail && "$(field threshold)" == 20 ]] || ok=1
expect threshold $ok "$(summary)"

run_check revcomp_hint "$WORKDIR/simple.fastq.gz" "$SIMPLE" 10 "$WORKDIR/wl16_rc.txt"
ok=0
[[ "$(field verdict)" == fail && "$(field pct_revcomp)" == 100.0 ]] || ok=1
[[ "$(field note)" == *"reverse complement"* ]] || ok=1
grep -q "reverse complement" "$LOG" || ok=1
expect revcomp_hint $ok "$(summary)"

run_check whitelist_width "$WORKDIR/simple.fastq.gz" "$SIMPLE" 10 "$WORKDIR/wl14.txt"
ok=0
[[ "$(field verdict)" == fail ]] || ok=1
[[ "$(field note)" == *"14-nt barcodes"* && "$(field note)" == *"16 nt"* ]] || ok=1
expect whitelist_width $ok "$(summary)"

# --------------------------------------------------------------------------
# Several barcode segments
# --------------------------------------------------------------------------

log_step "several segments"

run_check complex_pass "$WORKDIR/complex.fastq.gz" "$COMPLEX" 10 \
    "$WORKDIR/seg1.txt" "$WORKDIR/seg2.txt" "$WORKDIR/seg3.txt"
ok=0
[[ "$(field verdict)" == pass && "$(field pct_match)" == 100.0 ]] || ok=1
[[ "$(field pct_per_segment)" == "100.0,100.0,100.0" && "$(field layout_length)" == 40 ]] || ok=1
expect complex_pass $ok "$(summary)"

run_check complex_segment "$WORKDIR/complex.fastq.gz" "$COMPLEX" 10 \
    "$WORKDIR/seg1.txt" "$WORKDIR/segX.txt" "$WORKDIR/seg3.txt"
ok=0
[[ "$(field verdict)" == fail ]] || ok=1
IFS=, read -r p1 p2 p3 <<<"$(field pct_per_segment)"
[[ "$p1" == 100.0 && "$p3" == 100.0 ]] || ok=1
awk -v p="$p2" 'BEGIN { exit !(p < 5) }' || ok=1
grep -q "only segment(s) 2 fail" "$LOG" || ok=1
expect complex_segment $ok "$(summary)"

# --------------------------------------------------------------------------
# Nothing to check
# --------------------------------------------------------------------------

log_step "skipped"

run_check no_whitelist "$WORKDIR/simple.fastq.gz" "$SIMPLE" 10
ok=0
[[ $RC -eq 0 && "$(field verdict)" == skipped && "$(field note)" == *"no barcode whitelist"* ]] || ok=1
expect no_whitelist $ok "$(summary)"

run_check anchored_layout "$WORKDIR/complex.fastq.gz" \
    "--soloType CB_UMI_Complex --soloCBposition 0_0_0_7 2_-8_2_-1 --soloUMIposition 3_1_3_8" 10 \
    "$WORKDIR/seg1.txt" "$WORKDIR/seg2.txt"
ok=0
[[ $RC -eq 0 && "$(field verdict)" == skipped && "$(field note)" == *"2_-8_2_-1"* ]] || ok=1
expect anchored_layout $ok "$(summary)"

run_check whitelist_count "$WORKDIR/complex.fastq.gz" "$COMPLEX" 10 "$WORKDIR/seg1.txt"
ok=0
[[ "$(field verdict)" == skipped && "$(field note)" == "1 whitelist file(s) for 3 barcode segment(s)" ]] || ok=1
expect whitelist_count $ok "$(summary)"

run_check few_reads "$WORKDIR/few.fastq.gz" "$SIMPLE" 10 "$WORKDIR/wl16.txt"
ok=0
[[ "$(field verdict)" == skipped && "$(field note)" == *"only 200 reads"* ]] || ok=1
run_check few_reads_allowed "$WORKDIR/few.fastq.gz" "$SIMPLE" 10 --min-reads 100 "$WORKDIR/wl16.txt"
[[ "$(field verdict)" == pass ]] || ok=1
expect few_reads $ok "$(summary)"

# --------------------------------------------------------------------------
# Read length and input handling
# --------------------------------------------------------------------------

log_step "read length and input"

run_check short_read "$WORKDIR/short.fastq.gz" "$SIMPLE" 10 "$WORKDIR/wl16.txt"
ok=0
[[ "$(field verdict)" == fail && "$(field read_length)" == 26 ]] || ok=1
[[ "$(field note)" == *"26 nt"* && "$(field note)" == *"28 nt"* ]] || ok=1
grep -q "BARCODE CHECK FAILED" "$LOG" || ok=1
expect short_read $ok "$(summary)"

run_check long_read "$WORKDIR/long.fastq.gz" "$SIMPLE" 10 "$WORKDIR/wl16.txt"
ok=0
[[ "$(field verdict)" == pass && "$(field read_length)" == 30 ]] || ok=1
[[ "$(field note)" == *"longer than the 28 nt"* ]] || ok=1
expect long_read $ok "$(summary)"

run_check plain_input "$WORKDIR/simple.fastq" "$SIMPLE" 10 "$WORKDIR/wl16.txt.gz"
ok=0
[[ "$(field verdict)" == pass && "$(field pct_match)" == 100.0 ]] || ok=1
expect plain_input $ok "$(summary)"

# 1,500 reads on the whitelist, then 5,000 that are not: only the first are matched
run_check first_reads_only "$WORKDIR/head.fastq.gz" "$SIMPLE" 10 --reads 1500 "$WORKDIR/wl16.txt"
ok=0
[[ $RC -eq 0 && "$(field verdict)" == pass && "$(field reads)" == 1500 && "$(field pct_match)" == 100.0 ]] || ok=1
expect first_reads_only $ok "$(summary)"

run_check missing_file "$WORKDIR/none.fastq.gz" "$SIMPLE" 10 "$WORKDIR/wl16.txt"
[[ $RC -eq 2 && ! -e "$TSV" ]] && r=0 || r=1; expect missing_file $r "rc=$RC"

finish_check
