#!/usr/bin/env bash
#
# A read's round-1 barcode is matched at the smallest edit (Levenshtein) distance, 
# at most max_edit_dist, over A/C/G/T only. Between 8-nt strings, distance 2 includes 
# one deletion plus one insertion (a shifted barcode), and a barcode with N never matches. 
# The matches are looked up in a table built once from the whitelist; a scan of the whole 
# whitelist per read is kept as the reference:
#
#   - every barcode at 0-3 substitutions or with a deletion + insertion from a whitelist
#     with near-identical entries (ties), with Ns, lowercase and other characters;
#   - fixed examples of rules (shifted barcode, N, ties);
#   - well  notation: 'A1:C6' block, 'A1-B6' range in plate order, 'C4';
#   - end to end on a synthetic FASTQ pair: the same reads written, the same stats.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh parse_demux [-- OPTIONS]
       tests/checks/parse_demux.sh [OPTIONS]

Check parsebio_custom_demux.py's barcode matching (lookup table vs a whitelist scan,
 edit-distance rules) and its well notation.

Options:
  --python PATH   Python interpreter (default: \$BCA_PYTHON, else python3).
  --keep          Keep the generated fixture and outputs for inspection.
  -h, --help      Show this message.

Cases:
  table_matches   every barcode, ties, indels, Ns and odd characters get the scan's answer
  parse_rules     fixed examples: shifted barcode at distance 2, N never matches, ties kept
  well_ranges     well notation: blocks, plate-order ranges, single wells
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

def shift(s):
    """One deletion plus one insertion: at most 2 edits from s, same length."""
    i = rng.randrange(len(s)); short = s[:i] + s[i + 1:]
    j = rng.randrange(len(short) + 1)
    return short[:j] + rng.choice(B) + short[j:]

# A whitelist with near-identical entries (ties at distance 1 and 2), and a rotated copy
# of one entry (edit distance 2, Hamming distance up to 8).
wl = ["".join(rng.choice(B) for _ in range(8)) for _ in range(90)]
wl += [mutate(wl[0], 1), mutate(wl[1], 2), mutate(wl[2], 2), mutate(wl[2], 3), wl[4][1:] + wl[4][0]]
wl = list(dict.fromkeys(wl))

queries = set()
for w in wl:
    for n in range(4):
        queries.add(mutate(w, n))
    queries.add(mutate(w, 1, "ACGTN"))
    queries.add(shift(w)); queries.add(mutate(shift(w), 1))
for _ in range(2000):
    w = rng.choice(wl)
    queries.add(rng.choice([mutate(w, rng.randint(0, 3), "ACGTN"), shift(w), mutate(shift(w), 1)]))
queries |= {"".join(rng.choice(B) for _ in range(8)) for _ in range(1000)}
queries |= {"ACGT.ACG", wl[5].lower(), wl[0][:7] + "X", "NNNNNNNN", wl[0][:7] + "N", wl[3][:4], ""}

scan = {q: demux.best_matches_all(q, wl, max_ed=2) for q in queries}
mismatched = []
for ed in (0, 1, 2):
    t = demux.build_match_table(wl, max_ed=ed)
    for q in queries:
        want = scan[q] if scan[q][1] is not None and scan[q][1] <= ed else ([], None)
        if demux.best_matches(q, t) != want:
            mismatched.append((ed, q))
table = demux.build_match_table(wl, max_ed=2)
ties = sum(1 for d, seqs in table.values() if len(seqs) > 1)
indel_only = sum(1 for q, (best, d) in scan.items()
                 if best and min(sum(x != y for x, y in zip(q.upper(), b)) for b in best) > d)
with_n = [q for q in queries if "N" in q.upper()]
check("table_matches", not mismatched and ties > 0 and indel_only > 0 and all(not scan[q][0] for q in with_n),
      f"{len(queries)} barcodes x max_ed 0-2: {ties} tied table entries, {indel_only} matched only through "
      f"an indel, {len(with_n)} with N unmatched; {len(mismatched)} differ from the scan"
      + (f", e.g. {mismatched[:3]}" if mismatched else ""))

# Rules on fixed barcodes
w1, w2, w3 = "ACGTACGT", "TTGGCCAA", "TTGGCAAA"
t = demux.build_match_table([w1, w2, w3], max_ed=2)
cases = [
    ("exact",                      w1,         ([w1], 0)),
    ("one substitution",           "ACGTACGA", ([w1], 1)),
    ("base lost, next shifted in", "ACGACGTG", ([w1], 2)),   # T at 3 deleted, G appended
    ("rotated by one",             "CGTACGTA", ([w1], 2)),   # Hamming 8, edit distance 2
    ("three substitutions",        "AGGTTCGA", ([], None)),
    ("N in barcode",               "ACGTACGN", ([], None)),
    ("tie at distance 1",          "TTGGCGAA", ([w2, w3], 1)),
    ("lowercase",                  "acgtacgt", ([w1], 0)),
]
bad = [(n, q, demux.best_matches(q, t), want) for n, q, want in cases if demux.best_matches(q, t) != want]
check("parse_rules", not bad, f"{len(cases)} fixed examples" + (f"; wrong: {bad}" if bad else " as Parse's tables"))

# Well notation
row = lambda r, a, b: [f"{r}{c}" for c in range(a, b + 1)]
wcases = [
    ("A1:C6", row("A", 1, 6) + row("B", 1, 6) + row("C", 1, 6)),
    ("A1-B6", row("A", 1, 12) + row("B", 1, 6)),
    ("A7-B6", row("A", 7, 12) + row("B", 1, 6)),
    ("A1-A3", row("A", 1, 3)),
    ("A1:D12", [f"{r}{c}" for r in "ABCD" for c in range(1, 13)]),
    ("C4", ["C4"]),
    ("A1-A6,B1:D3,C4", row("A", 1, 6) + row("B", 1, 3) + row("C", 1, 4) + row("D", 1, 3)),
    ("a1:b2,A2", ["A1", "A2", "B1", "B2"]),
    ("B3:A1", []),
]
wbad = [(s, demux.parse_well_range(s), want) for s, want in wcases if demux.parse_well_range(s) != want]
for s in ("I1", "A13", "A0", "A1-", "XYZ"):
    try:
        demux.parse_well_range(s); wbad.append((s, "no error", "ValueError"))
    except ValueError:
        pass
check("well_ranges", not wbad, f"{len(wcases)} specifications and 5 malformed wells"
      + (f"; wrong: {wbad}" if wbad else " as Parse expands them"))

# End to end: the same reads, in the same order, and the same stats as the scan.
wells = [f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)]
wl_csv = os.path.join(work, "wl.csv")
with open(wl_csv, "w") as fh:
    fh.write("bci,sequence,uid,well,stype\n")
    for i, (s, w) in enumerate(zip(wl, wells)):
        fh.write(f"{i + 1},{s},u{i},{w},R\n")
fq1, fq2 = os.path.join(work, "R1.fastq.gz"), os.path.join(work, "R2.fastq.gz")
with gzip.open(fq1, "wt") as f1, gzip.open(fq2, "wt") as f2:
    for i in range(6000):
        w = rng.choice(wl)
        bc = shift(w) if i % 7 == 0 else mutate(w, rng.choice([0, 0, 1, 2, 3]), "ACGTN")
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
table_matches = demux.best_matches
demux.best_matches = lambda bc, t: demux.best_matches_all(bc, wl, max_ed=2)
try:
    scan_run = run("scan")
finally:
    demux.best_matches = table_matches
stats = table_run[0]
check("end_to_end", table_run == scan_run and stats["matched"] > 0 and stats["too_short"] > 0
      and stats["no_group_match"] > 0 and stats["n_in_bc"] > 0,
      f"{stats['total']} pairs, {stats['matched']} written, {stats['no_group_match']} other wells, "
      f"{stats['unmatched']} unmatched ({stats['n_in_bc']} with N), {stats['too_short']} too short; "
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
