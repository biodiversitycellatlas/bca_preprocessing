#!/usr/bin/env bash
# description: Check that MERGE_FASTQS concatenates gzipped FASTQs without recompressing them, and keeps every read in order
#
# MERGE_FASTQS (modules/local/custom/manipulate/merge_files) joins the files of one sample per
# read. A run of gzip members is itself a valid gzip file, so gzipped inputs are appended byte
# for byte and only plain inputs are compressed. The failures guarded against:
#   - a gzipped input decompressed and recompressed after all, which is what made the step slow;
#   - reads dropped, duplicated or reordered when gzipped and plain inputs are mixed;
#   - compression detected from the file name instead of the content;
#   - index reads (I1/I2) merged into the wrong file, or index outputs where there are none;
#   - a merged file that a reader stops in after the first gzip member.
#
# The module's own script runs, through a minimal Nextflow project on the local executor.
# Needs nextflow and pigz; the multi_member_read case also needs a python3.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

MODULE="$PROJECT_ROOT/modules/local/custom/manipulate/merge_files/main.nf"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh merge_fastqs [-- OPTIONS]
       tests/checks/merge_fastqs.sh [OPTIONS]

Run MERGE_FASTQS on small FASTQs and check what it writes.

Options:
  --python PATH   Python interpreter for multi_member_read (default: \$BCA_PYTHON,
                  else python3, else python).
  --keep          Keep the generated FASTQs and the Nextflow work directory.
  -h, --help      Show this message.

Cases:
  run                the module runs on every sample
  all_gz             three gzipped inputs: the output is their bytes, concatenated
  mixed              gzipped, plain, gzipped: the reads in order, gzip parts unchanged
  plain              two plain inputs: a gzipped output holding both
  single_gz          one gzipped input: symlinked
  mislabelled        compression follows the content, not the extension
  indices            I1 and I2 over two lanes go to their own outputs
  no_indices         no index outputs without index inputs
  multi_member_read  python's gzip reads every read of a multi-member output
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --python)  PYTHON="${2:?--python needs a value}"; shift 2 ;;
        --keep)    KEEP=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

# --------------------------------------------------------------------------
# Pre-flight
# --------------------------------------------------------------------------

log_header "merge_fastqs"

for tool in nextflow pigz; do
    if ! have_cmd "$tool"; then
        record SKIP "run" "$tool not found"
        finish_check
        exit 0
    fi
done

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"
    elif have_cmd python; then PYTHON="python"
    fi
fi

WORKDIR="$BCA_TEST_LOGDIR/merge_fastqs"
rm -rf "$WORKDIR"
mkdir -p "$WORKDIR/in"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

IN="$WORKDIR/in"

# --------------------------------------------------------------------------
# Fixture: small FASTQs whose read names say where each read came from
# --------------------------------------------------------------------------

# fastq NAME N LEN: N reads of LEN bases, named after NAME so the order is visible
fastq() {
    local name="$1" n="$2" len="$3" i bases qual
    bases="$(printf 'ACGT%.0s' $(seq 1 "$len") | cut -c1-"$len")"
    qual="$(printf 'F%.0s' $(seq 1 "$len"))"
    for ((i = 1; i <= n; i++)); do
        printf '@%s_%d\n%s\n+\n%s\n' "$name" "$i" "$bases" "$qual"
    done
}

# plain and gzipped versions of each part; the gzipped ones at different levels, so a
# recompressed output cannot match their bytes by chance
for part in A B C; do
    fastq "cdna_$part" 5 30 >"$IN/cdna_$part.fq"
    fastq "bc_$part" 5 16 >"$IN/bc_$part.fq"
done
for part in L1 L2; do
    for read in I1 I2; do
        fastq "${read}_$part" 4 8 >"$IN/${read}_$part.fq"
    done
done
level=1
for f in "$IN"/*.fq; do
    gzip -n -"$level" -c "$f" >"$f.gz"
    level=$(( level % 9 + 1 ))
done

# Sample all_gz: three gzipped files per read, no index reads
# Sample mixed: cDNA gzipped, plain, gzipped; BC two plain files (the 'plain' case)
# Sample single: one gzipped file per read
# Sample mislab: gzip named .fastq and plain named .fastq.gz; BC one gzip named .fastq
# Sample idx: two lanes, with I1 and I2
mkdir -p "$IN/mislab"
cp "$IN/cdna_A.fq.gz" "$IN/mislab/cdna_gz_named_plain.fastq"
cp "$IN/cdna_B.fq" "$IN/mislab/cdna_plain_named_gz.fastq.gz"
cp "$IN/bc_A.fq.gz" "$IN/mislab/bc_gz_named_plain.fastq"
mkdir -p "$IN/idx"
for part in L1 L2; do
    cp "$IN/I1_$part.fq.gz" "$IN/idx/S_${part}_I1.fastq.gz"
    cp "$IN/I2_$part.fq.gz" "$IN/idx/S_${part}_I2.fastq.gz"
done
printf 'samplesheet\n' >"$IN/input.csv"

# --------------------------------------------------------------------------
# Run the module
# --------------------------------------------------------------------------

NF="$WORKDIR/nf"
mkdir -p "$NF"
cat >"$NF/main.nf" <<NF
include { MERGE_FASTQS } from '$MODULE'

def f(names) {
    return names.collect { n -> file("$IN/\${n}") }
}

workflow {
    def sheet = file("$IN/input.csv")
    channel.of(
        [[id: 'all_gz'], f(['cdna_A.fq.gz', 'cdna_B.fq.gz', 'cdna_C.fq.gz']), f(['bc_A.fq.gz', 'bc_B.fq.gz', 'bc_C.fq.gz']), [], sheet],
        [[id: 'mixed'],  f(['cdna_A.fq.gz', 'cdna_B.fq', 'cdna_C.fq.gz']),    f(['bc_A.fq', 'bc_B.fq']),                     [], sheet],
        [[id: 'single'], f(['cdna_A.fq.gz']),                                 f(['bc_A.fq.gz']),                             [], sheet],
        [[id: 'mislab'], f(['mislab/cdna_gz_named_plain.fastq', 'mislab/cdna_plain_named_gz.fastq.gz']), f(['mislab/bc_gz_named_plain.fastq']), [], sheet],
        [[id: 'idx'],    f(['cdna_A.fq.gz', 'cdna_B.fq.gz']),                 f(['bc_A.fq.gz', 'bc_B.fq.gz']),
            f(['idx/S_L1_I1.fastq.gz', 'idx/S_L2_I1.fastq.gz', 'idx/S_L1_I2.fastq.gz', 'idx/S_L2_I2.fastq.gz']), sheet]
    ) | MERGE_FASTQS
    MERGE_FASTQS.out.merged_files.subscribe { meta, cdna, bc, ind, _sheet ->
        println "OUT\t\${meta.id}\t\${cdna}\t\${bc}\t\${[ind].flatten().join(',')}"
    }
}
NF

RUN_LOG="$BCA_TEST_LOGDIR/merge_fastqs_run.log"
if ! (cd "$NF" && run_logged "$RUN_LOG" nextflow -q run main.nf); then
    record FAIL "run" "nextflow run failed, see merge_fastqs_run.log"
    tail_log "$RUN_LOG"
    finish_check
    exit 1
fi
if [[ "$(grep -c '^OUT' "$RUN_LOG")" -ne 5 ]]; then
    record FAIL "run" "expected 5 merged samples, see merge_fastqs_run.log"
    finish_check
    exit 1
fi
record PASS "run" "5 samples merged"

# out SAMPLE FIELD: the output path of one sample (FIELD 3 = cDNA, 4 = BC, 5 = indices)
out() { awk -F'\t' -v s="$1" -v c="$2" '$1 == "OUT" && $2 == s { print $c }' "$RUN_LOG"; }

is_gzip() { [[ "$(od -An -N2 -tx1 "$1" | tr -d ' \n')" == "1f8b" ]]; }

# check NAME DETAIL CONDITION...: PASS when the command succeeds
check() {
    local name="$1" detail="$2"; shift 2
    if "$@"; then record PASS "$name" "$detail"; else record FAIL "$name" "$detail: not so"; fi
}

# --------------------------------------------------------------------------
# Cases
# --------------------------------------------------------------------------

all_gz() {
    local c b
    c="$(out all_gz 3)"; b="$(out all_gz 4)"
    cmp -s "$c" <(cat "$IN"/cdna_{A,B,C}.fq.gz) \
        && cmp -s "$b" <(cat "$IN"/bc_{A,B,C}.fq.gz) \
        && cmp -s <(gzip -dc "$c") <(cat "$IN"/cdna_{A,B,C}.fq) \
        && cmp -s <(gzip -dc "$b") <(cat "$IN"/bc_{A,B,C}.fq)
}
check "all_gz" "the inputs' bytes concatenated, and every read in order" all_gz

mixed() {
    local c sa sc
    c="$(out mixed 3)"
    sa="$(stat -c %s "$IN/cdna_A.fq.gz")"; sc="$(stat -c %s "$IN/cdna_C.fq.gz")"
    cmp -s <(gzip -dc "$c") <(cat "$IN"/cdna_{A,B,C}.fq) \
        && cmp -s <(head -c "$sa" "$c") "$IN/cdna_A.fq.gz" \
        && cmp -s <(tail -c "$sc" "$c") "$IN/cdna_C.fq.gz"
}
check "mixed" "reads in input order, the gzipped parts unchanged" mixed

plain() {
    local b
    b="$(out mixed 4)"
    is_gzip "$b" && cmp -s <(gzip -dc "$b") <(cat "$IN"/bc_{A,B}.fq)
}
check "plain" "plain inputs compressed into one gzipped output" plain

single_gz() {
    local c b
    c="$(out single 3)"; b="$(out single 4)"
    [[ -L "$c" && -L "$b" ]] && cmp -s "$c" "$IN/cdna_A.fq.gz" && cmp -s "$b" "$IN/bc_A.fq.gz"
}
check "single_gz" "a lone gzipped input is symlinked" single_gz

mislabelled() {
    local c b sa
    c="$(out mislab 3)"; b="$(out mislab 4)"
    sa="$(stat -c %s "$IN/cdna_A.fq.gz")"
    cmp -s <(gzip -dc "$c") <(cat "$IN"/cdna_{A,B}.fq) \
        && cmp -s <(head -c "$sa" "$c") "$IN/cdna_A.fq.gz" \
        && [[ -L "$b" ]] && cmp -s <(gzip -dc "$b") "$IN/bc_A.fq"
}
check "mislabelled" "gzip named .fastq copied, plain named .fastq.gz compressed" mislabelled

indices() {
    local ind i1 i2
    ind="$(out idx 5)"
    i1="$(tr ',' '\n' <<<"$ind" | grep '_merged_I1.fastq.gz$' || true)"
    i2="$(tr ',' '\n' <<<"$ind" | grep '_merged_I2.fastq.gz$' || true)"
    [[ -n "$i1" && -n "$i2" ]] \
        && cmp -s "$i1" <(cat "$IN"/idx/S_L{1,2}_I1.fastq.gz) \
        && cmp -s "$i2" <(cat "$IN"/idx/S_L{1,2}_I2.fastq.gz)
}
check "indices" "I1 and I2 each concatenated over both lanes" indices

no_indices() {
    local s
    for s in all_gz mixed single mislab; do
        [[ -z "$(out "$s" 5)" ]] || return 1
    done
}
check "no_indices" "no index outputs for the samples without index reads" no_indices

if [[ -z "$PYTHON" ]]; then
    record SKIP "multi_member_read" "no python found"
else
    got="$("$PYTHON" -c '
import gzip, sys
with gzip.open(sys.argv[1], "rt") as fh:
    print(sum(1 for _ in fh) // 4)
' "$(out all_gz 3)" 2>&1 || true)"
    if [[ "$got" == "15" ]]; then
        record PASS "multi_member_read" "15 reads over 3 gzip members"
    else
        record FAIL "multi_member_read" "expected 15 reads, got '$got'"
    fi
fi

finish_check
