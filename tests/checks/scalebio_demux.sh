#!/usr/bin/env bash
# description: ScaleBio demultiplexing (bin/scalebio_demux.py): barcode correction, sample split, trimming, synthetic barcode read
#
# bin/scalebio_demux.py corrects each barcode against its list (Hamming distance),
# drops ambiguous ones, trims the cDNA and splits the reads by RT (and PCR) well. Checked
# on synthetic reads with random barcodes (tests/lib/make_scalebio_fixture.py), whose
# truth.tsv holds every read's expected outcome:
#
#   - per layout (QuantumScale, v1.1) and i5 orientation: the read categories, and every
#     written read's sample, synthetic barcode read and trimmed length;
#   - the match table's rules: exact, corrected, ambiguous, N, too far;
#   - well notation: column-major plate ranges, ids holding dashes, malformed specs;
#   - barcode-list checks in 'prepare' and the padding of short barcodes;
#   - conf/seqtech_parameters.config reads the synthetic read where the script writes it;
#   - read files of different lengths are an error; 'gather' sums the chunk metrics.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh scalebio_demux [-- OPTIONS]
       tests/checks/scalebio_demux.sh [OPTIONS]

Check scalebio_demux.py on synthetic ScaleBio-layout reads.

Options:
  --python PATH   Python interpreter (default: \$BCA_PYTHON, else python3).
  --keep          Keep the generated fixture and outputs for inspection.
  -h, --help      Show this message.

Cases:
  quantum           QuantumScale reads: categories and every written read as expected
  quantum_i5_rc     the same with the i5 index reverse-complemented (detected from the reads)
  v1_1              v1.1 reads (splint-anchored barcodes, 8/9-nt ligation barcodes)
  v1_1_i5_rc        the same, orientation given rather than detected
  match_rules       exact / corrected / ambiguous / N / too far
  well_ranges       plate ranges run column by column; ids with dashes; bad specs rejected
  prepare_checks    duplicate sequences, too-long barcodes, wrong list count rejected; padding
  config_layout     seqtech_parameters' STARsolo and alevin-fry positions match the script
  uneven_reads      read files of different lengths make the run fail
  gather            a sample's metrics are the sum of its chunks
  decode            QuantumScale cell barcodes decode back to their sample's wells
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

log_header "scalebio_demux"

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"; elif have_cmd python; then PYTHON="python"; fi
fi
if [[ -z "$PYTHON" ]]; then
    record SKIP "tools" "a Python is needed"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/scalebio_demux"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

OUT="$WORKDIR/results.tsv"
if "$PYTHON" - "$PROJECT_ROOT" "$WORKDIR" "$PYTHON" >"$OUT" 2>"$WORKDIR/stderr.log" <<'PY'
import csv, gzip, json, os, re, shutil, subprocess, sys
from collections import Counter
root, work, python = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, os.path.join(root, "bin"))
import scalebio_demux as sb

SCRIPT = os.path.join(root, "bin", "scalebio_demux.py")
FIXTURE = os.path.join(root, "tests", "lib", "make_scalebio_fixture.py")

def check(name, ok, detail):
    print(("PASS" if ok else "FAIL") + "\t" + name + "\t" + detail.replace("\t", " ").replace("\n", " "))

def run(*args, ok=True):
    p = subprocess.run([python, *map(str, args)], capture_output=True, text=True)
    if ok and p.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, args[:2]))} failed: {p.stderr[-400:]}")
    return p

def fixture(name, layout, revcomp, reads=1400, seed=11):
    d = os.path.join(work, name)
    run(FIXTURE, "--layout", layout, "--outdir", d, "--reads", reads, "--seed", seed,
        *(["--i5-revcomp"] if revcomp else []))
    with open(os.path.join(d, "whitelists.txt")) as fh:
        lists = fh.read().split()
    run(SCRIPT, "prepare", "--layout", layout, "--whitelists", *lists, "--outdir", os.path.join(d, "prep"))
    return d

def demux(d, layout, chunk="lib1_01_L001", extra=(), ok=True, cdna=None):
    rd = lambda r: os.path.join(d, f"lib_S1_L001_{r}_001.fastq.gz")
    cd, bc = ("R1", "R2") if layout == "quantum" else ("R2", "R1")
    args = [SCRIPT, "demux", "--barcodes", os.path.join(d, "prep", "barcodes.json"),
            "--samples", os.path.join(d, "samples.tsv"), "--chunk", chunk,
            "--cdna", cdna or rd(cd), "--bc", rd(bc), "--i2", rd("I2"), "--outdir", os.path.join(d, "out")]
    if layout == "quantum":
        args += ["--i1", rd("I1")]
    return run(*args, *extra, ok=ok)

def verify(name, layout, revcomp, extra=()):
    try:
        d = fixture(name, layout, revcomp)
        demux(d, layout, extra=extra)
    except RuntimeError as exc:
        check(name, False, str(exc)); return
    truth = {r["read"]: r for r in csv.DictReader(open(os.path.join(d, "truth.tsv")), delimiter="\t")}
    want = Counter(r["category"] for r in truth.values())
    m = json.load(open(os.path.join(d, "out", "lib1_01_L001.metrics.json")))
    got = {k: v for k, v in m["read_categories"].items() if v}
    bad, seen = [], 0
    for s in ("sampleA", "sampleB", "sampleC"):
        with gzip.open(os.path.join(d, "out", f"{s}.lib1_01_L001_BC_UMI.fastq.gz"), "rt") as a, \
             gzip.open(os.path.join(d, "out", f"{s}.lib1_01_L001_cDNA.fastq.gz"), "rt") as b:
            la, lb = a.read().split("\n"), b.read().split("\n")
        for i in range(0, len(la) - 1, 4):
            seen += 1
            t = truth[la[i][1:]]
            if (t["sample"], t["bc_read"], int(t["cdna_len"])) != (s, la[i + 1], len(lb[i + 1])) \
                    or la[i] != lb[i] or len(lb[i + 1]) != len(lb[i + 3]):
                bad.append(la[i][1:])
    ok = got == dict(want) and seen == want["Pass"] and not bad and m["i5_revcomp"] == revcomp
    check(name, ok, f"{m['reads']} reads, {got}" + ("" if got == dict(want) else f" (expected {dict(want)})")
          + f"; {seen} written, {len(bad)} wrong; i5 reverse-complemented: {m['i5_revcomp']}")

verify("quantum", "quantum", False)
verify("quantum_i5_rc", "quantum", True)
verify("v1_1", "v1.1", False)
verify("v1_1_i5_rc", "v1.1", True, extra=("--i5-revcomp", "true"))

# Match rules: two entries 2 apart make their midpoints ambiguous at 1 substitution
seg = sb.Segment({"name": "x", "entries": [
    {"id": "a", "sequences": ["AAAAAAAA"], "canonical": "AAAAAAAA"},
    {"id": "b", "sequences": ["AAAAAACC"], "canonical": "AAAAAACC"},
    {"id": "c", "sequences": ["GGGGTTTT", "CCCCTTTT"], "canonical": "GGGGTTTT"}]}, 1)
cases = [("AAAAAAAA", 0), ("AAAATAAA", 0), ("AAAAAAAC", None), ("AAAAAACC", 1), ("AAAAAACG", 1),
         ("AAANAAAA", 0), ("AANNAAAA", None), ("TTAAAAAA", None), ("CCCCTTTA", 2), ("GGGG", None)]
wrong = [(o, e, seg.match(o)) for o, e in cases if seg.match(o) != e]
seg2 = sb.Segment({"name": "x", "entries": [{"id": "a", "sequences": ["AAAAAAAA"], "canonical": "AAAAAAAA"}]}, 2)
wrong += [("TTAAAAAA@2", 0, seg2.match("TTAAAAAA"))] if seg2.match("TTAAAAAA") != 0 else []
check("match_rules", not wrong, f"{len(cases) + 1} observed barcodes" + (f"; wrong: {wrong}" if wrong else " as expected"))

# Well notation
ids = [f"{c}{r}" for c in range(1, 13) for r in "ABCDEFGH"]
pcr = ["QSR-1", "QSR-2", "QSR-3", "QSR-10", "QSR-P"]
wcases = [
    (("1A-2C", ids), ids[:11]),
    (("2H-3B", ids), ["2H", "3A", "3B"]),
    (("1A;3B;5A-5C", ids), ["1A", "3B", "5A", "5B", "5C"]),
    (("", ids), ids),
    (("QSR-2-QSR-10", pcr), ["QSR-2", "QSR-3", "QSR-10"]),
    (("QSR-1;QSR-P", pcr), ["QSR-1", "QSR-P"]),
    (("1A 1A;1B", ids), ["1A", "1B"]),
]
wbad = [(a, sb.expand_wells(*a), e) for a, e in wcases if sb.expand_wells(*a) != e]
for s in ["13A", "2C-1A", "1A-", "X1", "QSR-3-QSR-1"]:
    try:
        sb.expand_wells(s, pcr if s.startswith("QSR") else ids); wbad.append((s, "accepted"))
    except ValueError:
        pass
check("well_ranges", not wbad, f"{len(wcases)} specifications, 5 malformed" + (f"; wrong: {wbad}" if wbad else " handled"))

# prepare: list checks and padding
pd = os.path.join(work, "prep_checks"); os.makedirs(pd, exist_ok=True)
def lst(name, lines):
    p = os.path.join(pd, name)
    with open(p, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return p
lig = lst("lig.tsv", ["L1\tACGTACGT", "L2\tTTGGCCAAG", "L3\tCCAATTGG"])
rt = lst("rt.tsv", [f"{w}\t{'ACGT'[i % 4]}{'ACGT'[i // 4 % 4]}{'ACGT'[i // 16 % 4]}AAAAAAA" for i, w in enumerate(ids[:8])])
pcr_ok = lst("pcr.tsv", ["P1\tACGTACGTAC", "P2\tTTTTGGGGCC"])
problems = []
p = run(SCRIPT, "prepare", "--layout", "v1.1", "--whitelists", lig, rt, pcr_ok, "--outdir", pd, ok=False)
if p.returncode != 0:
    problems.append(f"valid lists rejected: {p.stderr.strip()[-200:]}")
else:
    canon = open(os.path.join(pd, "whitelist_1_lig.txt")).read().split()
    if canon != ["ACGTACGTA", "TTGGCCAAG", "CCAATTGGA"]:
        problems.append(f"ligation barcodes padded to {canon}")
for label, lists in [
    ("duplicate sequence", [lig, rt, lst("dup.tsv", ["P1\tACGTACGTAC", "P2\tACGTACGTAC"])]),
    ("too long", [lig, rt, lst("long.tsv", ["P1\tACGTACGTACG"])]),
    ("not a sequence", [lig, rt, lst("bad.tsv", ["P1\tACGTXCGTAC"])]),
    ("two lists for three segments", [lig, rt]),
]:
    if run(SCRIPT, "prepare", "--layout", "v1.1", "--whitelists", *lists, "--outdir", pd, ok=False).returncode == 0:
        problems.append(f"{label} accepted")
check("prepare_checks", not problems, "valid lists read, ligation barcodes padded to 9 nt, 4 bad inputs rejected"
      if not problems else "; ".join(problems))

# conf/seqtech_parameters.config reads the synthetic read where the script writes it
conf = open(os.path.join(root, "conf", "seqtech_parameters.config")).read()
cproblems = []
for proto, layout in [("scalebio_quantum_rna", "quantum"), ("scalebio_rna_v1_1", "v1.1")]:
    block = re.search(r"'" + re.escape(proto) + r"'\s*\{(.*?)\n\s*\}", conf, re.S)
    if not block:
        cproblems.append(f"{proto}: no block"); continue
    block = block.group(1)
    lay = sb.LAYOUTS[layout]
    pos, start = [], 0
    for s in lay["segments"]:
        pos.append(f"0_{start}_0_{start + lay['width'][s] - 1}"); start += lay["width"][s]
    umi = f"0_{start}_0_{start + lay['umi'][2] - 1}"
    # With the per-segment whitelists the segments are read one by one (CB_UMI_Complex);
    # without, the whole corrected barcode is one CB (CB_UMI_Simple, no whitelist)
    with_wl = re.search(r"scalebio_mapping_whitelist\s*=\s*(true|false)", block)
    if not with_wl:
        cproblems.append(f"{proto}: no scalebio_mapping_whitelist"); continue
    if with_wl.group(1) == "true":
        want_star = f"--soloType CB_UMI_Complex --soloCBposition {' '.join(pos)} --soloUMIposition {umi}"
    else:
        start = sb.cb_length(layout)
        want_star = f"--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen {start} --soloUMIstart {start + 1} --soloUMIlen {lay['umi'][2]}"
        if start > 31:
            cproblems.append(f"{proto}: a {start}-nt cell barcode is longer than STARsolo's 31")
    if want_star not in block:
        cproblems.append(f"{proto}: STARsolo layout is not '{want_star}'")
    if f'alevin_bc_geometry =    "1[1-{start}]"' not in block or \
       f'alevin_umi_geometry =   "1[{start + 1}-{start + lay["umi"][2]}]"' not in block:
        cproblems.append(f"{proto}: alevin-fry geometry is not 1[1-{start}] / 1[{start + 1}-{start + lay['umi'][2]}]")
    if f'scalebio_layout =       "{layout}"' not in block:
        cproblems.append(f"{proto}: scalebio_layout is not '{layout}'")
check("config_layout", not cproblems, "both protocols read the synthetic read as written" if not cproblems else "; ".join(cproblems))

# Read files of different lengths
d = os.path.join(work, "quantum")
short = os.path.join(work, "short_R1.fastq.gz")
with gzip.open(os.path.join(d, "lib_S1_L001_R1_001.fastq.gz"), "rt") as fi, gzip.open(short, "wt") as fo:
    fo.writelines(line for n, line in enumerate(fi) if n < 4 * 1000)
shutil.rmtree(os.path.join(d, "out"), ignore_errors=True)
p = demux(d, "quantum", ok=False, cdna=short)
check("uneven_reads", p.returncode != 0 and "different numbers of reads" in p.stderr,
      "a cDNA file 400 reads short " + ("fails the run" if p.returncode != 0 else "was ACCEPTED"))

# gather: two chunks of one library
try:
    shutil.rmtree(os.path.join(d, "out"), ignore_errors=True)
    demux(d, "quantum", chunk="lib1_01_L001")
    demux(d, "quantum", chunk="lib1_02_L002")
    g = os.path.join(work, "gathered.json")
    run(SCRIPT, "gather", "--sample", "sampleB", "--metrics",
        os.path.join(d, "out", "lib1_01_L001.metrics.json"), os.path.join(d, "out", "lib1_02_L002.metrics.json"), "--out", g)
    one = json.load(open(os.path.join(d, "out", "lib1_01_L001.metrics.json")))
    out = json.load(open(g))
    ok = (out["sample_reads"] == 2 * one["samples"]["sampleB"] and out["library_reads"] == 2 * one["reads"]
          and out["library_read_categories"]["Pass"] == 2 * one["read_categories"]["Pass"] and len(out["chunks"]) == 2)
    check("gather", ok, f"sampleB: {out['sample_reads']} reads over 2 chunks of {one['samples']['sampleB']}")
except RuntimeError as exc:
    check("gather", False, str(exc))

# decode: the QuantumScale cell barcode gives back each read's barcode ids
import random
rng = random.Random(3)
sizes = [96, 96, 96, 96, 13]
trip = [[rng.randrange(s) for s in sizes] for _ in range(2000)] + [[s - 1 for s in sizes], [0] * 5]
round_trip = all(sb.decode_cell(sb.encode_cell(t, sizes, 16), sizes) == t for t in trip)
cells = os.path.join(work, "cells.txt")
wells = {"sampleA": ([f"{c}{r}" for c in range(1, 13) for r in "ABCDEFGH"][:11], None),
         "sampleB": (["3A", "3B"] + [f"5{r}" for r in "ABCDEFGH"], [f"P-{i}" for i in range(1, 7)]),
         "sampleC": (["3A", "3B"], [f"P-{i}" for i in range(7, 13)])}
dbad, n_cells = [], 0
for s, (rts, pcrs) in wells.items():
    with gzip.open(os.path.join(d, "out", f"{s}.lib1_01_L001_BC_UMI.fastq.gz"), "rt") as fh:
        cbs = sorted({line[:16] for n, line in enumerate(fh) if n % 4 == 1})
    with open(cells, "w") as fh:
        fh.write("\n".join(cbs) + "\n")
    p = run(SCRIPT, "decode", "--barcodes", os.path.join(d, "prep", "barcodes.json"), "--cells", cells)
    rows = [l.split("\t") for l in p.stdout.strip().split("\n")[1:]]
    n_cells += len(rows)
    dbad += [r for r in rows if r[4] not in rts or (pcrs and r[5] not in pcrs) or "?" in r]
check("decode", round_trip and not dbad and n_cells > 0,
      f"2002 random cells round-trip {'exactly' if round_trip else 'WRONG'}; "
      f"{n_cells} written cells decode to their sample's wells" + (f"; {len(dbad)} do not: {dbad[:2]}" if dbad else ""))
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
