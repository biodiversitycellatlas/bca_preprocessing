#!/usr/bin/env bash
# description: Check that PARSEBIO_CUSTOM_DEMUX's barcode lookup table gives exactly the per-read whitelist scan's matches
#
# parsebio_custom_demux.py assigns each read to wells by the whitelist barcodes at the
# smallest Hamming distance (<= max_edit_dist). Scanning the whole whitelist for every
# read made it the slowest step of a Parse run, so the matches are now looked up in a
# table built once from the whitelist. The table must answer exactly as the scan did,
# ties included, so the scan is kept as the reference:
#
#   - every barcode at 0-3 mismatches from a whitelist with near-identical entries
#     (ties at distance 1 and 2), with Ns, and with characters outside ACGTN, which
#     the table does not cover and leaves to the scan;
#   - end to end on a synthetic FASTQ pair: the same reads written, the same stats.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh parse_demux [-- OPTIONS]
       tests/checks/parse_demux.sh [OPTIONS]

Compare parsebio_custom_demux.py's barcode lookup table with the whitelist scan it
replaces, per barcode and end to end on a synthetic FASTQ pair.

Options:
  --python PATH   Python interpreter (default: \$BCA_PYTHON, else python3).
  --keep          Keep the generated fixture and outputs for inspection.
  -h, --help      Show this message.

Cases:
  table_matches   every barcode, ties, Ns and odd characters get the scan's answer
  end_to_end      the same reads are written, with the same stats
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

log_header "parse_demux"

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"; elif have_cmd python; then PYTHON="python"; fi
fi
if [[ -z "$PYTHON" ]]; then
    record SKIP "tools" "a Python is needed"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/parse_demux"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

OUT="$WORKDIR/results.tsv"
if "$PYTHON" - "$PROJECT_ROOT/bin" "$WORKDIR" >"$OUT" 2>"$WORKDIR/stderr.log" <<'PY'
import gzip, itertools, os, random, sys
sys.path.insert(0, sys.argv[1])
import parsebio_custom_demux as demux

work = sys.argv[2]
rng = random.Random(7)
B = "ACGT"

def check(name, ok, detail):
    print(("PASS" if ok else "FAIL") + "\t" + name + "\t" + detail)

def mutate(s, n, alphabet=B):
    s = list(s)
    for p in rng.sample(range(len(s)), n):
        s[p] = rng.choice([c for c in alphabet if c != s[p]])
    return "".join(s)

# A whitelist with near-identical entries, so ties at distance 1 and 2 occur.
wl = ["".join(rng.choice(B) for _ in range(8)) for _ in range(90)]
wl += [mutate(wl[0], 1), mutate(wl[1], 2), mutate(wl[2], 2), mutate(wl[2], 3)]
wl = list(dict.fromkeys(wl))
table = demux.build_match_table(wl, max_ed=2)

queries = set()
for w in wl:
    for n in range(4):
        queries.add(mutate(w, n))
        queries.add(mutate(w, n, "ACGTN"))
for _ in range(20000):
    w = rng.choice(wl)
    queries.add(mutate(w, rng.randint(0, 3), "ACGTN"))
queries |= {"".join(rng.choice(B) for _ in range(8)) for _ in range(5000)}
queries |= {"ACGT.ACG", "acgtacgt", wl[0][:7] + "X", "NNNNNNNN", wl[3][:4]}

mismatched = []
for ed in (0, 1, 2):
    t = demux.build_match_table(wl, max_ed=ed)
    for q in queries:
        if demux.best_matches(q, t, wl, max_ed=ed) != demux.best_matches_all(q, wl, max_ed=ed):
            mismatched.append((ed, q))
ties = sum(1 for d, seqs in table.values() if len(seqs) > 1)
check("table_matches", not mismatched and ties > 0,
      f"{len(queries)} barcodes x max_ed 0-2, {ties} tied table entries; "
      f"{len(mismatched)} differ from the scan" + (f", e.g. {mismatched[:3]}" if mismatched else ""))

# End to end: the same reads, in the same order, and the same stats as the scan.
wells = [f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)]
wl_csv = os.path.join(work, "wl.csv")
with open(wl_csv, "w") as fh:
    fh.write("bci,sequence,uid,well,stype\n")
    for i, (s, w) in enumerate(zip(wl, wells)):
        fh.write(f"{i + 1},{s},u{i},{w},R\n")
fq1, fq2 = os.path.join(work, "R1.fastq.gz"), os.path.join(work, "R2.fastq.gz")
with gzip.open(fq1, "wt") as f1, gzip.open(fq2, "wt") as f2:
    for i in range(20000):
        bc = mutate(rng.choice(wl), rng.choice([0, 0, 1, 2, 3]), "ACGTN")
        length = 86 if i % 50 else 55   # some reads too short to hold the barcode
        s2 = ("".join(rng.choice(B) for _ in range(50)) + bc +
              "".join(rng.choice(B) for _ in range(28)))[:length]
        s1 = "".join(rng.choice(B) for _ in range(60))
        f1.write(f"@r{i}\n{s1}\n+\n{'I' * 60}\n")
        f2.write(f"@r{i}\n{s2}\n+\n{'I' * len(s2)}\n")

wl_full = demux.load_whitelist_csv(wl_csv)
group = set(demux.parse_well_range("A1:D12"))

def run(tag):
    o1, n1 = demux.open_out(work, tag, "G", "R1")
    o2, n2 = demux.open_out(work, tag, "G", "R2")
    try:
        stats = demux.process_fastq(fq1, fq2, wl_full, group, 50, 58, o1, o2, max_ed=2)
    finally:
        o1.close(); o2.close()
    with gzip.open(n1, "rt") as a, gzip.open(n2, "rt") as b:
        return stats, a.read(), b.read()

table_run = run("table")
scan_matches = demux.best_matches
demux.best_matches = lambda bc, t, w, max_ed=2: demux.best_matches_all(bc, w, max_ed=max_ed)
try:
    scan_run = run("scan")
finally:
    demux.best_matches = scan_matches
stats = table_run[0]
check("end_to_end", table_run == scan_run and stats["matched"] > 0 and stats["too_short"] > 0
      and stats["no_group_match"] > 0,
      f"{stats['total']} pairs, {stats['matched']} written, {stats['no_group_match']} other wells, "
      f"{stats['unmatched']} unmatched, {stats['too_short']} too short; "
      + ("identical to the scan" if table_run == scan_run else "DIFFERS from the scan"))
PY
then
    while IFS=$'\t' read -r status name detail; do
        record "$status" "$name" "$detail"
    done <"$OUT"
    [[ -s "$OUT" ]] || record FAIL "python" "no results produced"
else
    record FAIL "python" "$(tail -3 "$WORKDIR/stderr.log" 2>/dev/null | tr '\n' ' ')"
fi

finish_check
