#!/usr/bin/env bash
# description: Check the SPLiT-seq barcode lists, RT-well demultiplexing and round-1 merge on synthetic v1 reads
#
# Academic SPLiT-seq (Rosenberg et al. 2018) holds UMI 0-10, round 3 10-18, round 2 48-56
# and round 1 (RT) 86-94 in read 2 (0-based). Each RT well holds oligo-dT primer i (1-48)
# and random-hexamer primer i+48. This check runs the scripts behind SPLITSEQ_WELLS,
# SPLITSEQ_DEMUX and SPLITSEQ_COLLAPSE_BC1 with the positions the 'splitseq_v1' block gives:
#
#   - splitseq_wells.py derives the well CSV from the 96 barcodes (i and i+48 share well
#     i, A1-D12, T/R) and reads lists bare, delimited, headed or gzipped, as one or three;
#   - it refuses lists that are not 96 barcodes, rows without one barcode and duplicates;
#   - the demultiplexing keeps a sample's wells, oligo-dT and random-hexamer reads alike;
#   - the merge rewrites random-hexamer barcodes to their well's oligo-dT barcode, and
#     stops when most barcode reads are too short to hold round 1;
#   - parsebio_custom_demux.py's chemistry templates put round 1 at 86 (v1), 78 (v2), 50 (v3);
#   - the 'splitseq_v1' STARsolo, alevin-fry and demultiplexing positions agree.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh splitseq [-- OPTIONS]
       tests/checks/splitseq.sh [OPTIONS]

Run the SPLiT-seq pre-processing scripts on synthetic barcode lists and reads.

Options:
  --python PATH   Python interpreter (default: \$BCA_PYTHON, else python3).
  --keep          Keep the generated fixture and outputs for inspection.
  -h, --help      Show this message.

Cases:
  wells_derived   96 barcodes -> well CSV (i and i+48 share well i, A1-D12, T/R) and three plain lists
  list_formats    bare, gzipped, '#'-headed table, header row and three lists give the same lists
  list_refused    95 barcodes, a row with two barcodes, a duplicate and two files are refused
  wells_given     a given well CSV is copied as it is and the lists still written
  demux_v1        rt 'A1-A6' keeps the oligo-dT and random-hexamer reads of A1-A6 (0-1 edits), nothing else
  collapse_v1     random-hexamer round-1 barcodes become their well's oligo-dT barcode; whitelist = 48 oligo-dT
  too_short       a 66 nt barcode read stops the merge with a message
  chemistry       --chemistry v1/v2/v3 gives round 1 at 86-94 / 78-86 / 50-58
  seqtech_v1      'splitseq_v1' STARsolo positions = alevin-fry geometry = splitseq_bc1_start/end
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

log_header "splitseq"

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"; elif have_cmd python; then PYTHON="python"; fi
fi
if [[ -z "$PYTHON" ]]; then
    record SKIP "tools" "a Python is needed"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/splitseq"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

OUT="$WORKDIR/results.tsv"
if "$PYTHON" - "$PROJECT_ROOT" "$WORKDIR" >"$OUT" 2>"$WORKDIR/stderr.log" <<'PY'
import csv, gzip, os, random, re, subprocess, sys
root, work = sys.argv[1], sys.argv[2]
bindir = os.path.join(root, "bin")
sys.path.insert(0, bindir)
rng = random.Random(7)
B = "ACGT"
env = dict(os.environ, PYTHONPATH=bindir)

def check(name, ok, detail):
    print(("PASS" if ok else "FAIL") + "\t" + name + "\t" + detail)

def run(script, *args):
    return subprocess.run([sys.executable, os.path.join(bindir, script), *map(str, args)],
                          capture_output=True, text=True, env=env)

def last_err(proc):
    lines = (proc.stderr or proc.stdout).strip().splitlines()
    return lines[-1] if lines else f"exit {proc.returncode}"

def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]

def mutate(s, n):
    s = list(s)
    for p in rng.sample(range(len(s)), n):
        s[p] = rng.choice([c for c in B if c != s[p]])
    return "".join(s)

def write(path, text, gz=False):
    with (gzip.open(path, "wt") if gz else open(path, "w")) as fh:
        fh.write(text)
    return path

def lines(path):
    with open(path) as fh:
        return fh.read().split()

# 96 barcodes at least 3 edits apart, so a read 1 edit off is closest to its own only
# (distance-2 matching and ties are covered by parse_demux and parse_collapse)
bcs = []
while len(bcs) < 96:
    s = "".join(rng.choice(B) for _ in range(8))
    if all(lev(s, b) >= 3 for b in bcs):
        bcs.append(s)

# --- wells_derived ------------------------------------------------------------------
lst = write(os.path.join(work, "bc_8nt_v1.csv"), "".join(f"{i + 1},{s}\n" for i, s in enumerate(bcs)))
out = os.path.join(work, "wells")
p = run("splitseq_wells.py", lst, "--output", out)
if p.returncode:
    check("wells_derived", False, last_err(p))
else:
    with open(os.path.join(out, "splitseq_wells.csv")) as fh:
        rows = list(csv.DictReader(fh))
    plate = [f"{r}{c}" for r in "ABCD" for c in range(1, 13)]
    ok = (len(rows) == 96
          and all(r["sequence"] == bcs[i] and r["bci"] == str(i + 1) for i, r in enumerate(rows))
          and all(rows[i]["well"] == rows[i + 48]["well"] == plate[i] for i in range(48))
          and all(r["stype"] == ("T" if i < 48 else "R") for i, r in enumerate(rows))
          and all(lines(os.path.join(out, f"splitseq_bc{k}.txt")) == bcs for k in (1, 2, 3)))
    check("wells_derived", ok, f"{len(rows)} rows; barcodes i and i+48 in well i (A1..D12), T/R; three lists of 96")

# --- list_formats -------------------------------------------------------------------
variants = {
    "bare":    write(os.path.join(work, "bare.txt"), "\n".join(bcs) + "\n"),
    "gzipped": write(os.path.join(work, "bare.txt.gz"), "\n".join(bcs) + "\n", gz=True),
    "table":   write(os.path.join(work, "table.tsv"), "#WellPosition\tName\tSequence\n"
                     + "".join(f"{'ABCDEFGH'[i // 12]}{i % 12 + 1}\tRound1_{i + 1:02d}\t{s.lower()}\n" for i, s in enumerate(bcs))),
    "header":  write(os.path.join(work, "header.csv"), "index,sequence\n" + "".join(f"{i + 1},{s}\n" for i, s in enumerate(bcs))),
}
bad = []
for name, path in variants.items():
    o = os.path.join(work, f"fmt_{name}")
    p = run("splitseq_wells.py", path, "--output", o)
    if p.returncode or any(lines(os.path.join(o, f"splitseq_bc{k}.txt")) != bcs for k in (1, 2, 3)):
        bad.append(f"{name}: {last_err(p) if p.returncode else 'lists differ'}")
r2 = write(os.path.join(work, "r2.txt"), "\n".join(bcs[:24]) + "\n")
r3 = write(os.path.join(work, "r3.txt"), "\n".join(bcs[24:40]) + "\n")
o = os.path.join(work, "fmt_three")
p = run("splitseq_wells.py", variants["bare"], r2, r3, "--output", o)
if p.returncode or [lines(os.path.join(o, f"splitseq_bc{k}.txt")) for k in (1, 2, 3)] != [bcs, bcs[:24], bcs[24:40]]:
    bad.append(f"three lists: {last_err(p) if p.returncode else 'lists differ'}")
check("list_formats", not bad, "; ".join(bad) if bad else "bare, gzipped, '#'-headed table, header row, three lists")

# --- list_refused -------------------------------------------------------------------
cases = {
    "95 barcodes": ([write(os.path.join(work, "short.txt"), "\n".join(bcs[:95]) + "\n")], "95 round-1 barcodes"),
    "two in a row": ([write(os.path.join(work, "two.csv"), "".join(f"{s},{bcs[0]}\n" for s in bcs))], "found 2"),
    "duplicate": ([write(os.path.join(work, "dup.txt"), "\n".join(bcs[:95] + [bcs[3]]) + "\n")], "listed twice"),
    "two files": ([variants["bare"], variants["bare"]], "gives 2 files"),
}
bad = []
for name, (args, msg) in cases.items():
    p = run("splitseq_wells.py", *args, "--output", os.path.join(work, "refused"))
    if p.returncode == 0 or msg not in (p.stderr + p.stdout):
        bad.append(f"{name}: exit {p.returncode}, {last_err(p)}")
check("list_refused", not bad, "; ".join(bad) if bad else ", ".join(cases))

# --- wells_given --------------------------------------------------------------------
given = write(os.path.join(work, "my_wells.csv"), "bci,sequence,uid,well,stype\n"
              + "".join(f"{i + 1},{s},u{i},H{i % 12 + 1},T\n" for i, s in enumerate(bcs[:12])))
o = os.path.join(work, "given")
p = run("splitseq_wells.py", r2, "--wells", given, "--output", o)
ok = (p.returncode == 0 and open(os.path.join(o, "splitseq_wells.csv")).read() == open(given).read()
      and lines(os.path.join(o, "splitseq_bc1.txt")) == bcs[:24])
check("wells_given", ok, "copied as given; 24-barcode list kept" if ok else last_err(p))

# --- synthetic v1 reads ---------------------------------------------------------------
L1, L2 = "GTGGCCGATGTTTCGCATCGGCGTACGACT", "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"
def r2_read(bc1, length=110):
    umi = "".join(rng.choice(B) for _ in range(10))
    s = umi + rng.choice(bcs) + L1 + rng.choice(bcs) + L2 + bc1 + "T" * 30
    assert s[86:94] == bc1
    return s[:length]

wells_csv = os.path.join(work, "wells", "splitseq_wells.csv")
group = {f"A{c}" for c in range(1, 7)}
recs, keep = [], []
for _ in range(2000):
    i = rng.randrange(96)
    recs.append(r2_read(mutate(bcs[i], rng.choice([0, 0, 1]))))
    keep.append((i % 48) < 6)                           # wells A1-A6, oligo-dT or random hexamer
for _ in range(30):                                     # far from every barcode: dropped
    while True:
        s = "".join(rng.choice(B) for _ in range(8))
        if min(lev(s, b) for b in bcs) > 2:
            break
    recs.append(r2_read(s)); keep.append(False)
fq1 = os.path.join(work, "S1_R1.fastq.gz"); fq2 = os.path.join(work, "S1_R2.fastq.gz")
with gzip.open(fq1, "wt") as f1, gzip.open(fq2, "wt") as f2:
    for n, r in enumerate(recs):
        f1.write(f"@r{n} 1:N:0:1\n{'ACGT' * 16}\n+\n{'F' * 64}\n")
        f2.write(f"@r{n} 2:N:0:1\n{r}\n+\n{'F' * len(r)}\n")

# --- demux_v1 -------------------------------------------------------------------------
o = os.path.join(work, "demux")
p = run("parsebio_custom_demux.py", "--sample_id", "S1", "--fq1", fq1, "--fq2", fq2, "--whitelist", wells_csv,
        "--group", "S1", "A1-A6", "--output", o, "--barcode_start", 86, "--barcode_end", 94, "--max_edit_dist", 2)
if p.returncode:
    check("demux_v1", False, last_err(p))
else:
    with gzip.open(os.path.join(o, "S1_group_S1_R2.fastq.gz"), "rt") as fh:
        got = [l.split()[0][1:] for k, l in enumerate(fh) if k % 4 == 0]
    want = [f"r{n}" for n, k in enumerate(keep) if k]
    check("demux_v1", got == want, f"{len(got)} of {len(recs)} pairs kept, {len(want)} expected (A1-A6, T and R, 0-1 edits)")

# --- collapse_v1 ----------------------------------------------------------------------
o = os.path.join(work, "collapse")
p = run("parsebio_collapse_bc1.py", "--csv", wells_csv, "--fq2", fq2, "--sample_id", "S1", "--output", o,
        "--barcode_start", 86, "--barcode_end", 94, "--max_edit_dist", 2,
        "--whitelist_name", "splitseq_bc1_T_whitelist.txt", "--max_too_short", 0.5)
if p.returncode:
    check("collapse_v1", False, last_err(p))
else:
    with gzip.open(os.path.join(o, "S1_R2.bc1collapsed.fastq.gz"), "rt") as fh:
        seqs = [l.strip() for k, l in enumerate(fh) if k % 4 == 1]
    wrong = 0
    for orig, new in zip(recs[:2000], seqs):
        best = min(range(96), key=lambda j: lev(orig[86:94], bcs[j]))
        wrong += new[86:94] != bcs[best % 48] or new[:86] != orig[:86] or new[94:] != orig[94:]
    unchanged = all(a == b for a, b in zip(recs[2000:], seqs[2000:]))
    wl = lines(os.path.join(o, "splitseq_bc1_T_whitelist.txt"))
    ok = len(seqs) == len(recs) and not wrong and unchanged and wl == bcs[:48]
    check("collapse_v1", ok, f"{len(seqs)} reads, {wrong} wrong round-1 rewrites, far barcodes unchanged: {unchanged}; "
          f"whitelist {len(wl)} oligo-dT barcodes")

# --- too_short ------------------------------------------------------------------------
fq_short = os.path.join(work, "short_R2.fastq.gz")
with gzip.open(fq_short, "wt") as fh:
    for n in range(50):
        fh.write(f"@s{n}\n{r2_read(bcs[n], length=66)}\n+\n{'F' * 66}\n")
p = run("parsebio_collapse_bc1.py", "--csv", wells_csv, "--fq2", fq_short, "--sample_id", "S2",
        "--output", os.path.join(work, "short"), "--barcode_start", 86, "--barcode_end", 94, "--max_too_short", 0.5)
check("too_short", p.returncode != 0 and "shorter than 94 nt" in p.stderr,
      "stopped: " + last_err(p) if p.returncode else "exit 0 on 66 nt barcode reads")

# --- chemistry ------------------------------------------------------------------------
from parsebio_custom_demux import bc1_bounds_from_chem
got = {c: bc1_bounds_from_chem(c) for c in ("v1", "v2", "v3")}
want = {"v1": (86, 94), "v2": (78, 86), "v3": (50, 58)}
check("chemistry", got == want, f"{got}")

# --- seqtech_v1 -----------------------------------------------------------------------
conf = open(os.path.join(root, "conf", "seqtech_parameters.config")).read()
m = re.search(r"'splitseq_v1'\s*\{(.*?)\n\s*\}", conf, re.S)
block = m.group(1) if m else ""
def key(k):
    mm = re.search(rf"\b{k}\s*=\s*(\"[^\"]*\"|\d+)", block)
    return mm.group(1).strip('"') if mm else None
try:
    solo = key("star_soloTypestring")
    cb = [tuple(map(int, x.split("_")[1::2])) for x in re.search(r"--soloCBposition ((?:\S+ ?)+?)(?= --)", solo).group(1).split()]
    umi = tuple(map(int, re.search(r"--soloUMIposition (\S+)", solo).group(1).split("_")[1::2]))
    alevin_cb = [tuple(map(int, x.split("-"))) for x in re.search(r"1\[(.*)\]", key("alevin_bc_geometry")).group(1).split(",")]
    alevin_umi = tuple(map(int, re.search(r"1\[(.*)\]", key("alevin_umi_geometry")).group(1).split("-")))
    bc1 = (int(key("splitseq_bc1_start")), int(key("splitseq_bc1_end")))
    # STAR: 0-based inclusive; alevin: 1-based inclusive; splitseq_bc1: 0-based half-open
    ok = (cb == [(86, 93), (48, 55), (10, 17)] and umi == (0, 9)
          and [(s + 1, e + 1) for s, e in cb] == alevin_cb and (umi[0] + 1, umi[1] + 1) == alevin_umi
          and (cb[0][0], cb[0][1] + 1) == bc1 and key("cdna_read") == "R1")
    check("seqtech_v1", ok, f"STAR CB {cb} UMI {umi}; alevin {alevin_cb} {alevin_umi}; bc1 {bc1}")
except Exception as e:
    check("seqtech_v1", False, f"block not parsed: {e!r}")
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
