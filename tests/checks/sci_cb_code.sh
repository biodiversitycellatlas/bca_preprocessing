#!/usr/bin/env bash
# description: Check the sciRNAseq3 cell barcode code that alevin-fry maps (bin/sci_cb_code.sh)
#
# sciRNAseq3's 40-nt cell barcode is longer than salmon alevin takes, so ALEVIN_FRY rewrites
# read 1 as a 16-nt code of the cell's place in the sample's four barcode lists, and decodes
# the matrix's barcodes back to STARsolo's p7_p5_lig_rt form. This check generates lists and
# reads and tests the round trip and the edge cases; no mapper is needed.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

CODE="$PROJECT_ROOT/bin/sci_cb_code.sh"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh sci_cb_code [-- OPTIONS]
       tests/checks/sci_cb_code.sh [OPTIONS]

Validate bin/sci_cb_code.sh against generated barcode lists and reads, with each awk
found on the PATH (awk, mawk, gawk, busybox awk).

Options:
  --keep          Keep the generated fixtures for inspection.
  -h, --help      Show this message.

Cases (per awk):
  round_trip       every read decodes to its own p7_p5_lig_rt; UMI and qualities kept
  read_layout      the rewritten read is code (16) + UMI (8), its qualities as long
  permit_list      the permit list is exactly the distinct codes in the reads
  known_codes      first entries give AAAA..., last entries give the largest code
  distinct_cells   different cells never share a code
  duplicates       a repeated list entry counts once
  unknown_segment  a barcode on no list, and a short read, get a code of N's, counted
  too_many         lists whose combinations exceed the code stop with exit 3
  entry_length     a list entry of another length than its segment stops with exit 3
  decode_errors    a code of the wrong length, or beyond the lists, stops with exit 1
  compressed       gzipped lists and reads give the same output as plain ones
  process_substitution        read through <(...) as ALEVIN_FRY does: same output, status by wait $!
  process_substitution_error  a failing encoder gives an empty stream and its status by wait $!
  no_reads         an empty read file gives an empty permit list
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep)    KEEP=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

log_header "sci_cb_code"

WORKDIR="$BCA_TEST_LOGDIR/sci_cb_code"
rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

# make_list FILE N LENGTH SEED -- N distinct random barcodes
make_list() {
    awk -v n="$2" -v len="$3" -v seed="$4" 'BEGIN {
        srand(seed); split("A C G T", B, " ")
        while (k < n) {
            s = ""; for (i = 0; i < len; i++) s = s B[int(rand() * 4) + 1]
            if (!(s in seen)) { seen[s] = 1; print s; k++ }
        }
    }' >"$1"
}

make_list "$WORKDIR/p7.txt"  3 10 1
make_list "$WORKDIR/p5.txt"  4 10 2
make_list "$WORKDIR/lig10.txt" 6 10 3
make_list "$WORKDIR/lig9.txt"  2  9 4
make_list "$WORKDIR/rt.txt" 12 10 5
# Ligation list as sci-rocket writes it: the 10-nt barcodes, then the 9-nt ones padded with G
{ cat "$WORKDIR/lig10.txt"; sed 's/$/G/' "$WORKDIR/lig9.txt"; } >"$WORKDIR/lig.txt"
LISTS=("$WORKDIR/p7.txt" "$WORKDIR/p5.txt" "$WORKDIR/lig.txt" "$WORKDIR/rt.txt")
SEGS="1-10,11-20,21-30,31-40"

# make_reads N SEED [LISTS...] -- 48-nt barcode reads (4 x 10 nt + 8-nt UMI); prints the
# truth (barcode joined with '_' and UMI) to FILE.truth
make_reads() {
    local n="$1" seed="$2" out="$3"; shift 3
    awk -v n="$n" -v seed="$seed" -v out="$out" -v truth="$out.truth" '
        FNR == 1 { f++ }
        { wl[f, FNR] = $1; wn[f] = FNR }
        END {
            srand(seed); split("A C G T", B, " ")
            for (r = 1; r <= n; r++) {
                s = ""; t = ""
                for (k = 1; k <= f; k++) { b = wl[k, int(rand() * wn[k]) + 1]; s = s b; t = t (k > 1 ? "_" : "") b }
                u = ""; for (k = 0; k < 8; k++) u = u B[int(rand() * 4) + 1]
                q = ""; for (k = 0; k < 48; k++) q = q substr("FGHI#?", int(rand() * 6) + 1, 1)
                printf "@r%d extra\n%s%s\n+\n%s\n", r, s, u, q > out
                print t "\t" u "\t" substr(q, 41, 8) > truth
            }
        }' "$@"
}

make_reads 3000 11 "$WORKDIR/reads.fastq" "${LISTS[@]}"
gzip -c "$WORKDIR/reads.fastq" >"$WORKDIR/reads.fastq.gz"
for f in "${LISTS[@]}"; do gzip -c "$f" >"$f.gz"; done

# The awks to test the script with: each one found on the PATH, put first as 'awk'
AWKS=()
for a in awk mawk gawk "busybox awk"; do
    bin=${a%% *}
    have_cmd "$bin" || continue
    [[ "$a" == "busybox awk" ]] && { busybox awk 'BEGIN{}' 2>/dev/null || continue; }
    AWKS+=("$a")
done

expect() {  # expect NAME RESULT DETAIL
    if [[ "$2" -eq 0 ]]; then record PASS "$1"; else record FAIL "$1" "$3"; fi
}

for awk_name in "${AWKS[@]}"; do
    tag=$(slugify "$awk_name")
    W="$WORKDIR/$tag"; mkdir -p "$W/bin"
    if [[ "$awk_name" == "busybox awk" ]]; then
        printf '#!/bin/sh\nexec busybox awk "$@"\n' >"$W/bin/awk"
    else
        printf '#!/bin/sh\nexec %s "$@"\n' "$(command -v "$awk_name")" >"$W/bin/awk"
    fi
    chmod +x "$W/bin/awk"
    run() { PATH="$W/bin:$PATH" bash "$CODE" "$@"; }

    log_step "with $awk_name"

    RC=0
    run encode --lists "${LISTS[@]}" --segments "$SEGS" --umi 41-48 --digits 16 \
        --r1 "$WORKDIR/reads.fastq" --out "$W/enc.fastq" --permit "$W/permit.txt" 2>"$W/enc.log" || RC=$?
    awk 'NR % 4 == 2 { print substr($0, 1, 16) }' "$W/enc.fastq" >"$W/codes.txt"
    run decode --lists "${LISTS[@]}" --digits 16 <"$W/codes.txt" >"$W/decoded.txt" 2>"$W/dec.log" || RC=$?
    awk 'NR % 4 == 2 { print substr($0, 17) } NR % 4 == 0 { print substr($0, 17) }' "$W/enc.fastq" \
        | paste - - >"$W/umi.txt"
    paste "$W/decoded.txt" "$W/umi.txt" >"$W/got.txt"
    ok=0; [[ $RC -eq 0 ]] && cmp -s "$W/got.txt" "$WORKDIR/reads.fastq.truth" || ok=1
    expect "round_trip:$tag" $ok "rc=$RC; $(diff "$W/got.txt" "$WORKDIR/reads.fastq.truth" | head -3 | tr '\n' ' ')"

    ok=0
    awk 'NR % 4 == 2 && length($0) != 24 { bad++ } NR % 4 == 0 && length($0) != 24 { bad++ }
         NR % 4 == 1 && $0 !~ /^@r[0-9]+ extra$/ { bad++ } NR % 4 == 3 && $0 != "+" { bad++ }
         END { exit bad > 0 }' "$W/enc.fastq" || ok=1
    [[ $(( $(wc -l <"$W/enc.fastq") / 4 )) -eq 3000 ]] || ok=1
    expect "read_layout:$tag" $ok "$(sed -n 1,4p "$W/enc.fastq" | tr '\n' '|')"

    ok=0
    sort -u "$W/codes.txt" >"$W/codes_u.txt"
    sort "$W/permit.txt" | cmp -s - "$W/codes_u.txt" || ok=1
    [[ $(wc -l <"$W/permit.txt") -eq $(sort -u "$W/permit.txt" | wc -l) ]] || ok=1
    grep -q "3000 reads, $(wc -l <"$W/permit.txt") cell barcodes, 0 reads with a barcode on no list" "$W/enc.log" || ok=1
    expect "permit_list:$tag" $ok "$(cat "$W/enc.log")"

    # First entries: code 0. Last entries: 3*4*8*12 - 1 = 1151, in 16 base-4 digits
    first="$(head -1 "$WORKDIR/p7.txt")$(head -1 "$WORKDIR/p5.txt")$(head -1 "$WORKDIR/lig.txt")$(head -1 "$WORKDIR/rt.txt")"
    last="$(tail -1 "$WORKDIR/p7.txt")$(tail -1 "$WORKDIR/p5.txt")$(tail -1 "$WORKDIR/lig.txt")$(tail -1 "$WORKDIR/rt.txt")"
    printf '@a\n%sAAAAAAAA\n+\n%s\n@b\n%sCCCCCCCC\n+\n%s\n' "$first" "$(printf 'I%.0s' {1..48})" "$last" "$(printf 'I%.0s' {1..48})" >"$W/known.fastq"
    run encode --lists "${LISTS[@]}" --segments "$SEGS" --umi 41-48 --digits 16 \
        --r1 "$W/known.fastq" --out "$W/known_enc.fastq" --permit "$W/known_permit.txt" 2>/dev/null
    n=1151; want=""; bases=ACGT
    for _ in $(seq 16); do want="${bases:$((n % 4)):1}$want"; n=$((n / 4)); done
    got=$(awk 'NR % 4 == 2 { print substr($0, 1, 16) }' "$W/known_enc.fastq" | paste -sd' ')
    [[ "$got" == "AAAAAAAAAAAAAAAA $want" ]] && r=0 || r=1
    expect "known_codes:$tag" $r "got '$got', want 'AAAAAAAAAAAAAAAA $want'"

    ok=0
    a=$(cut -f1 "$WORKDIR/reads.fastq.truth" | sort -u | wc -l); b=$(sort -u "$W/codes.txt" | wc -l)
    [[ "$a" -eq "$b" ]] || ok=1
    [[ $(paste "$W/codes.txt" <(cut -f1 "$WORKDIR/reads.fastq.truth") | sort -u | wc -l) -eq "$a" ]] || ok=1
    expect "distinct_cells:$tag" $ok "$a cells, $b codes"

    { cat "$WORKDIR/p5.txt"; head -2 "$WORKDIR/p5.txt"; } >"$W/p5_dup.txt"
    out=$(run check --lists "$WORKDIR/p7.txt" "$W/p5_dup.txt" "$WORKDIR/lig.txt" "$WORKDIR/rt.txt" --segments "$SEGS" --digits 16 2>&1) && r=0 || r=1
    [[ "$out" == *"(3 x 4 x 8 x 12 = 1152 combinations)"* ]] || r=1
    run encode --lists "$WORKDIR/p7.txt" "$W/p5_dup.txt" "$WORKDIR/lig.txt" "$WORKDIR/rt.txt" --segments "$SEGS" \
        --umi 41-48 --digits 16 --r1 "$WORKDIR/reads.fastq" --out "$W/dup_enc.fastq" --permit "$W/dup_permit.txt" 2>/dev/null || r=1
    cmp -s "$W/dup_enc.fastq" "$W/enc.fastq" || r=1
    expect "duplicates:$tag" $r "$out"

    # A read whose RT is on no list, then one of 30 nt
    { printf '@x\n%sTTTTTTTTTTGGGGGGGG\n+\n%s\n' "${first:0:30}" "$(printf 'I%.0s' {1..48})"
      printf '@y\n%s\n+\n%s\n' "${first:0:30}" "$(printf 'I%.0s' {1..30})"
      sed -n 1,4p "$WORKDIR/reads.fastq"; } >"$W/bad.fastq"
    # The fake RT must not be on the list by chance
    grep -qx TTTTTTTTTT "$WORKDIR/rt.txt" && log_warn "TTTTTTTTTT is on the RT list; unknown_segment is not meaningful"
    RC=0
    run encode --lists "${LISTS[@]}" --segments "$SEGS" --umi 41-48 --digits 16 \
        --r1 "$W/bad.fastq" --out "$W/bad_enc.fastq" --permit "$W/bad_permit.txt" 2>"$W/bad.log" || RC=$?
    ok=0; [[ $RC -eq 0 ]] || ok=1
    [[ "$(sed -n 2p "$W/bad_enc.fastq")" == "NNNNNNNNNNNNNNNNGGGGGGGG" ]] || ok=1
    [[ "$(sed -n 6p "$W/bad_enc.fastq")" == "NNNNNNNNNNNNNNNNNNNNNNNN" ]] || ok=1
    [[ "$(sed -n 8p "$W/bad_enc.fastq")" == "FFFFFFFFFFFFFFFF########" ]] || ok=1
    [[ $(wc -l <"$W/bad_permit.txt") -eq 1 ]] || ok=1
    grep -q "3 reads, 1 cell barcodes, 2 reads with a barcode on no list" "$W/bad.log" || ok=1
    expect "unknown_segment:$tag" $ok "rc=$RC $(cat "$W/bad.log") $(sed -n '2p;6p;8p' "$W/bad_enc.fastq" | tr '\n' ' ')"

    # 3 x 4 x 8 x 12 = 1152 combinations do not fit 5 digits (1024)
    RC=0; out=$(run check --lists "${LISTS[@]}" --segments "$SEGS" --digits 5 2>&1) || RC=$?
    ok=0; [[ $RC -eq 3 && "$out" == *"1152 barcode combinations, more than a 5-digit code holds (1024)"* ]] || ok=1
    RC2=0; run encode --lists "${LISTS[@]}" --segments "$SEGS" --umi 41-48 --digits 5 \
        --r1 "$WORKDIR/reads.fastq" --out "$W/tm.fastq" --permit "$W/tm.txt" 2>/dev/null || RC2=$?
    [[ $RC2 -eq 3 ]] || ok=1
    expect "too_many:$tag" $ok "check rc=$RC encode rc=$RC2: $out"

    { cat "$WORKDIR/rt.txt"; echo ACGTACGT; } >"$W/rt_short.txt"
    RC=0; out=$(run check --lists "$WORKDIR/p7.txt" "$WORKDIR/p5.txt" "$WORKDIR/lig.txt" "$W/rt_short.txt" \
        --segments "$SEGS" --digits 16 2>&1) || RC=$?
    [[ $RC -eq 3 && "$out" == *"list 4 has a 8-nt entry"* ]] && r=0 || r=1
    expect "entry_length:$tag" $r "rc=$RC: $out"

    ok=0
    RC=0; echo ACGT | run decode --lists "${LISTS[@]}" --digits 16 >/dev/null 2>"$W/de1.log" || RC=$?
    [[ $RC -eq 1 ]] && grep -q "not a 16-digit code" "$W/de1.log" || ok=1
    RC=0; echo TTTTTTTTTTTTTTTT | run decode --lists "${LISTS[@]}" --digits 16 >/dev/null 2>"$W/de2.log" || RC=$?
    [[ $RC -eq 1 ]] && grep -q "beyond the barcode lists" "$W/de2.log" || ok=1
    expect "decode_errors:$tag" $ok "$(cat "$W/de1.log" "$W/de2.log" | tr '\n' ' ')"

    run encode --lists "${LISTS[@]/%/.gz}" --segments "$SEGS" --umi 41-48 --digits 16 \
        --r1 "$WORKDIR/reads.fastq.gz" --out "$W/gz_enc.fastq" --permit "$W/gz_permit.txt" 2>/dev/null \
        && cmp -s "$W/gz_enc.fastq" "$W/enc.fastq" && cmp -s "$W/gz_permit.txt" "$W/permit.txt" && r=0 || r=1
    expect "compressed:$tag" $r "gzipped input gave other output"

    # As ALEVIN_FRY reads it: a process substitution, then 'wait $!' for the encoder's status
    cat <(run encode --lists "${LISTS[@]}" --segments "$SEGS" --umi 41-48 --digits 16 \
        --r1 "$WORKDIR/reads.fastq.gz" --out - --permit "$W/ps_permit.txt" 2>/dev/null) >"$W/ps_enc.fastq"
    RC=0; wait $! || RC=$?
    [[ $RC -eq 0 ]] && cmp -s "$W/ps_enc.fastq" "$W/enc.fastq" && cmp -s "$W/ps_permit.txt" "$W/permit.txt" && r=0 || r=1
    expect "process_substitution:$tag" $r "rc=$RC, or its output differs"

    # An encoder that fails gives its reader an empty stream at once, and 'wait $!' its status
    cat <(run encode --lists "${LISTS[@]}" --segments "$SEGS" --umi 41-48 --digits 5 \
        --r1 "$WORKDIR/reads.fastq.gz" --out - --permit "$W/ps_bad.txt" 2>/dev/null) >"$W/ps_bad.fastq"
    RC=0; wait $! || RC=$?
    [[ $RC -eq 3 && ! -s "$W/ps_bad.fastq" ]] && r=0 || r=1
    expect "process_substitution_error:$tag" $r "rc=$RC (want 3), $(wc -c <"$W/ps_bad.fastq") bytes read"

    : >"$W/empty.fastq"
    RC=0; run encode --lists "${LISTS[@]}" --segments "$SEGS" --umi 41-48 --digits 16 \
        --r1 "$W/empty.fastq" --out "$W/empty_enc.fastq" --permit "$W/empty_permit.txt" 2>/dev/null || RC=$?
    [[ $RC -eq 0 && -e "$W/empty_permit.txt" && ! -s "$W/empty_permit.txt" && ! -s "$W/empty_enc.fastq" ]] && r=0 || r=1
    expect "no_reads:$tag" $r "rc=$RC"
done

[[ ${#AWKS[@]} -gt 0 ]] || record FAIL "awk" "no awk on the PATH"

finish_check
