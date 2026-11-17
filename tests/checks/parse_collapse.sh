#!/usr/bin/env bash
# description: Check that PARSEBIO_COLLAPSE_BC1 merges each Parse well's T/R/X round-1 barcodes into its T barcode
#
# Each Parse round-1 well holds an oligo-dT (T), a random-hexamer (R) and in some kits
# targeted (X) RT primers, each with its own barcode. parsebio_collapse_bc1.py rewrites
# bc1 in the barcode read to the well's T barcode, so the mappers count one cell per well
# combination, as split-pipe does:
#
#   - T, R and X barcodes of a well, exact, with 1-2 mismatches or shifted by a deletion
#     plus an insertion (edit distance 2, as Parse's tables), become the well's T;
#   - reads tied between wells, unmatched (bc1 with N included) and too short keep their barcode;
#   - nothing else in a record changes, and no record is lost;
#   - the T whitelist holds one barcode per well in CSV order (first row if no T row).

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh parse_collapse [-- OPTIONS]
       tests/checks/parse_collapse.sh [OPTIONS]

Run parsebio_collapse_bc1.py on a synthetic round-1 CSV and barcode read.

Options:
  --python PATH   Python interpreter (default: \$BCA_PYTHON, else python3).
  --keep          Keep the generated fixture and outputs for inspection.
  -h, --help      Show this message.

Cases:
  collapse_map    T/R/X barcodes (exact, 1-2 mismatches, shifted) become the well's T; ties, N, unmatched, short reads unchanged
  records_intact  same records, headers, qualities and bases outside bc1
  t_whitelist     one barcode per well in CSV order, the first row for a well without a T barcode
  stats           the log's counts add up to the reads
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

log_header "parse_collapse"

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"; elif have_cmd python; then PYTHON="python"; fi
fi
if [[ -z "$PYTHON" ]]; then
    record SKIP "tools" "a Python is needed"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/parse_collapse"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

OUT="$WORKDIR/results.tsv"
if "$PYTHON" - "$PROJECT_ROOT/bin" "$WORKDIR" >"$OUT" 2>"$WORKDIR/stderr.log" <<'PY'
import gzip, os, random, subprocess, sys
bindir, work = sys.argv[1], sys.argv[2]
rng = random.Random(11)
B = "ACGT"

def check(name, ok, detail):
    print(("PASS" if ok else "FAIL") + "\t" + name + "\t" + detail)

def ham(a, b):
    return sum(x != y for x, y in zip(a, b))

def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]

def shift(s):
    """One deletion plus one insertion: at most 2 edits from s, same length."""
    i = rng.randrange(len(s)); short = s[:i] + s[i + 1:]
    j = rng.randrange(len(short) + 1)
    return short[:j] + rng.choice(B) + short[j:]

def mutate(s, n):
    s = list(s)
    for p in rng.sample(range(len(s)), n):
        s[p] = rng.choice([c for c in B if c != s[p]])
    return "".join(s)

# Two barcodes of different wells 2 apart (positions 6-7): a read 1 off both is a tie
# between wells. All other round-1 barcodes are at least 5 edits from each other and from
# this pair, so reads at most 2 edits off stay closest to their own barcode.
tie_a = "".join(rng.choice(B) for _ in range(8))
tie_b = tie_a[:6] + "".join(rng.choice([c for c in B if c != x]) for x in tie_a[6:])
bcs = []
for _ in range(300000):                      # seed 11 needs about 100,000 draws
    s = "".join(rng.choice(B) for _ in range(8))
    if all(ham(s, b) >= 5 and lev(s, b) >= 5 for b in bcs + [tie_a, tie_b]):
        bcs.append(s)
        if len(bcs) == 28:
            break
assert len(bcs) == 28, f"only {len(bcs)} barcodes 5 edits apart; change the seed"
wells = [f"A{i}" for i in range(1, 13)]
rows = []                                    # (bci, seq, well, stype)
for i, w in enumerate(wells):
    rows.append((i + 1, bcs[i], w, "T"))
for i, w in enumerate(wells):
    rows.append((97 + i, bcs[12 + i], w, "R"))
rows += [(200, bcs[24], "A4", "X1"), (231, bcs[25], "A10", "X4")]
# B1 has no T row: its first row in bci order (bci 150) is its canonical barcode.
rows += [(160, bcs[27], "B1", "R"), (150, bcs[26], "B1", "X_ibh")]
rows += [(240, tie_a, "C1", "T"), (241, tie_b, "C2", "T")]
csv_path = os.path.join(work, "bc_data_test_R1.csv")
with open(csv_path, "w") as fh:
    fh.write("bci,sequence,uid,well,stype\n")
    fh.writelines(f"{b},{s},u{b},{w},{t}\n" for b, s, w, t in rows)
seq_well = {s: w for _, s, w, _ in rows}
canon = {}
for b, s, w, t in sorted(rows):
    if t == "T":
        canon.setdefault(w, s)
for b, s, w, t in sorted(rows):
    canon.setdefault(w, s)

# Barcode read: UMI(10) bc3(8) linker(12) bc2(8) linker(12) bc1(8) + tail.
def read_with(bc1, length=86):
    return ("".join(rng.choice(B) for _ in range(50)) + bc1 + "".join(rng.choice(B) for _ in range(28)))[:length]

recs, expect = [], []                       # expected bc1 after the rewrite (None = unchanged)
for _ in range(3000):
    b, s, w, t = rng.choice([r for r in rows if r[2] not in ("C1", "C2")])
    n = rng.choice([0, 0, 1, 2])
    recs.append(read_with(mutate(s, n))); expect.append(canon[w])
shifted = 0
for _ in range(300):                         # a base lost inside bc1, the next one shifted in
    b, s, w, t = rng.choice([r for r in rows if r[2] not in ("C1", "C2")])
    r = shift(s)
    shifted += ham(r, s) > 2
    recs.append(read_with(r)); expect.append(canon[w])
for _ in range(20):                          # an N in bc1 never matches, as in Parse's tables
    s = rng.choice(bcs)
    p = rng.randrange(8)
    recs.append(read_with(s[:p] + "N" + s[p + 1:])); expect.append(None)
tie_read = tie_a[:6] + "".join(c for c in tie_b[6:7]) + tie_a[7]   # 1 off tie_a, 1 off tie_b
assert ham(tie_read, tie_a) == 1 and ham(tie_read, tie_b) == 1, (tie_read, tie_a, tie_b)
for _ in range(40):
    recs.append(read_with(tie_read)); expect.append(None)
far = []
while len(far) < 40:
    s = "".join(rng.choice(B) for _ in range(8))
    if min(lev(s, b) for b in seq_well) > 2:
        far.append(s)
for s in far:
    recs.append(read_with(s)); expect.append(None)
for _ in range(20):
    recs.append(read_with(bcs[0], length=55)); expect.append(None)
order = list(range(len(recs))); rng.shuffle(order)
recs = [recs[i] for i in order]; expect = [expect[i] for i in order]
quals = ["".join(rng.choice("FI:#,") for _ in r) for r in recs]

fq2 = os.path.join(work, "S1_R2.fastq.gz")
with gzip.open(fq2, "wt") as fh:
    for i, (r, q) in enumerate(zip(recs, quals)):
        fh.write(f"@r{i} 2:N:0:1\n{r}\n+\n{q}\n")

outdir = os.path.join(work, "out")
env = dict(os.environ, PYTHONPATH=bindir)
proc = subprocess.run([sys.executable, os.path.join(bindir, "parsebio_collapse_bc1.py"), "--csv", csv_path,
                       "--fq2", fq2, "--sample_id", "S1", "--output", outdir],
                      capture_output=True, text=True, env=env)
if proc.returncode:
    check("run", False, proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else f"exit {proc.returncode}")
    sys.exit()

with gzip.open(os.path.join(outdir, "S1_R2.bc1collapsed.fastq.gz"), "rt") as fh:
    lines = fh.read().splitlines()
out = [lines[i:i + 4] for i in range(0, len(lines), 4)]

wrong = [(i, recs[i][50:58], o[1][50:58], e) for i, (o, e) in enumerate(zip(out, expect))
         if o[1][50:58] != (e if e is not None else recs[i][50:58])]
n_rw = sum(e is not None for e in expect)
check("collapse_map", len(out) == len(recs) and not wrong,
      f"{n_rw} T/R/X reads (0-2 mismatches, {shifted} shifted beyond 2 mismatches) to their well's T barcode, "
      f"120 left unchanged (ties, N, unmatched, short); {len(wrong)} wrong" + (f", e.g. {wrong[:3]}" if wrong else ""))

intact = len(out) == len(recs) and all(
    o[0] == f"@r{i} 2:N:0:1" and o[2] == "+" and o[3] == quals[i] and len(o[1]) == len(recs[i])
    and o[1][:50] == recs[i][:50] and o[1][58:] == recs[i][58:]
    for i, o in enumerate(out))
check("records_intact", intact, f"{len(out)} of {len(recs)} records; headers, qualities and bases outside bc1 unchanged")

with open(os.path.join(outdir, "parse_bc1_T_whitelist.txt")) as fh:
    wl = fh.read().split()
want = []
for _, s, w, _ in rows:                      # wells in CSV order of first appearance
    if canon[w] not in want:
        want.append(canon[w])
check("t_whitelist", wl == want,
      f"{len(wl)} barcodes for {len(canon)} wells; B1 (no T row) -> bci 150" if wl == want else f"got {wl}, want {want}")

stats = {}
with open(os.path.join(outdir, "S1_bc1_collapse.log")) as fh:
    for l in fh:
        k, _, v = l.partition(":")
        stats[k] = v.strip()
parts = sum(int(stats[k]) for k in ("rewritten", "ambiguous", "unmatched", "too_short"))
froms = sum(int(v) for k, v in stats.items() if k.startswith("from_"))
eds = sum(int(v) for k, v in stats.items() if k.startswith("ed"))
ok = (int(stats["total"]) == len(recs) == parts and froms == eds == int(stats["rewritten"]) == n_rw
      and stats["ambiguous"] == "40" and stats["unmatched"] == "60" and stats["too_short"] == "20")
check("stats", ok, f"total {stats['total']} = rewritten {stats['rewritten']} + ambiguous {stats['ambiguous']} "
      f"+ unmatched {stats['unmatched']} + too short {stats['too_short']}; by type/distance sum to rewritten")
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
