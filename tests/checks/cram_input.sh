#!/usr/bin/env bash
# description: CRAM input: CRAM_TO_FASTQ round-trips Illumina-style CRAMs; ScaleBio demultiplexes Ultima CRAMs from their tags
#
# Illumina-style CRAMs (both reads, index reads in the BC/QT tags) are converted to FASTQ by
# CRAM_TO_FASTQ (modules/local/tools/samtools/cram_to_fastq), run here through a minimal
# Nextflow project. Ultima CRAMs, whose ScaleBio barcodes the Ultima trimmer moved into
# tags, are read by 'scalebio_demux.py demux --cram'. Fixture: random-barcode ScaleBio reads
# from tests/lib/make_scalebio_fixture.py, whose truth.tsv gives every read's outcome.
#
#   - dual index: FASTQs -> samtools import -> CRAM -> CRAM_TO_FASTQ gives the same reads,
#     index reads included, named like bcl-convert output;
#   - single index: one index read comes back as I1;
#   - unpaired: a CRAM of single-end reads fails instead of losing them;
#   - Ultima, UMI in a tag / at the start of the read: the same categories and written reads
#     as the FASTQ truth;
#   - tag values: barcode ids, ids inside read-group names, sequences.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

MODULE="$PROJECT_ROOT/modules/local/tools/samtools/cram_to_fastq/main.nf"
PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh cram_input [-- OPTIONS]
       tests/checks/cram_input.sh [OPTIONS]

Check CRAM input: CRAM_TO_FASTQ (Illumina-style CRAMs) and the Ultima CRAM demultiplexing.

Options:
  --python PATH   Python interpreter (default: \$BCA_PYTHON, else python3).
  --keep          Keep the fixtures, CRAMs and the Nextflow work directory.
  -h, --help      Show this message.

Cases:
  illumina_dual     dual-index CRAM -> FASTQ: identical reads and index reads, bcl-convert names
  illumina_single   single-index CRAM -> FASTQ: the index read comes back as I1
  unpaired          a single-end CRAM fails the conversion
  ultima_tag        Ultima CRAM, UMI in a tag: categories and reads as the FASTQ truth
  ultima_seq        Ultima CRAM, UMI at the start of the read: the same
  tag_values        barcode ids, ids inside read-group names and sequences are matched

Needs samtools and a Python; the illumina_* and unpaired cases also need nextflow.
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

log_header "cram_input"

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"; elif have_cmd python; then PYTHON="python"; fi
fi
if [[ -z "$PYTHON" ]] || ! have_cmd samtools; then
    record SKIP "tools" "samtools and a Python are needed"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/cram_input"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

FIX="$WORKDIR/fixture"
"$PYTHON" "$PROJECT_ROOT/tests/lib/make_scalebio_fixture.py" --layout quantum --outdir "$FIX/q" --reads 1400 --seed 21
"$PYTHON" "$PROJECT_ROOT/tests/lib/make_scalebio_fixture.py" --layout v1.1 --outdir "$FIX/v" --reads 300 --seed 22
Q="$FIX/q/lib_S1_L001"
V="$FIX/v/lib_S1_L001"

# --------------------------------------------------------------------------
# Illumina-style CRAMs through CRAM_TO_FASTQ
# --------------------------------------------------------------------------

CR="$WORKDIR/crams"
mkdir -p "$CR"
samtools import -1 "${Q}_R1_001.fastq.gz" -2 "${Q}_R2_001.fastq.gz" --i1 "${Q}_I1_001.fastq.gz" --i2 "${Q}_I2_001.fastq.gz" \
    -O cram -o "$CR/dual_S1_L001.cram"
samtools import -1 "${V}_R1_001.fastq.gz" -2 "${V}_R2_001.fastq.gz" --i1 "${V}_I2_001.fastq.gz" \
    -O cram -o "$CR/single_S1_L001.cram"
samtools import -0 "${Q}_R1_001.fastq.gz" -O cram -o "$CR/unpaired.cram"

if have_cmd nextflow; then
    NF="$WORKDIR/nf"
    mkdir -p "$NF"
    cp -r "$PROJECT_ROOT/lib" "$NF/lib"
    cat >"$NF/main.nf" <<NF
include { CRAM_TO_FASTQ } from '$MODULE'

params.cram_index_format = 'auto'
params.dynamic_time = null

workflow {
    def crams = channel.of(
        [[id: 'dual'],     file("$CR/dual_S1_L001.cram")],
        [[id: 'single'],   file("$CR/single_S1_L001.cram")],
        [[id: 'unpaired'], file("$CR/unpaired.cram")]
    )
    CRAM_TO_FASTQ(crams, [])
    CRAM_TO_FASTQ.out.fastqs.subscribe { meta, r1, r2, ind ->
        println "OUT\t\${meta.id}\t\${r1}\t\${r2}\t\${[ind].flatten().join(',')}"
    }
}
NF
    cat >"$NF/nextflow.config" <<'CFG'
process.errorStrategy = 'ignore'
CFG
    RUN_LOG="$BCA_TEST_LOGDIR/cram_input_run.log"
    (cd "$NF" && run_logged "$RUN_LOG" nextflow -q run main.nf) || true
fi

OUT="$WORKDIR/results.tsv"
if "$PYTHON" - "$PROJECT_ROOT" "$WORKDIR" "$PYTHON" "${RUN_LOG:-}" >"$OUT" 2>"$WORKDIR/stderr.log" <<'PY'
import csv, gzip, json, os, subprocess, sys
from collections import Counter
root, work, python, run_log = sys.argv[1:5]
sys.path.insert(0, os.path.join(root, "bin"))
import scalebio_demux as sb
SCRIPT = os.path.join(root, "bin", "scalebio_demux.py")
FIX = os.path.join(work, "fixture")

def check(name, ok, detail):
    print(("PASS" if ok else "FAIL") + "\t" + name + "\t" + detail.replace("\t", " ").replace("\n", " "))

def skip(name, detail):
    print("SKIP\t" + name + "\t" + detail)

def reads(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().split("\n")
    return [(lines[i].split()[0].lstrip("@"), lines[i + 1], lines[i + 3]) for i in range(0, len(lines) - 1, 4)]

# ---- CRAM_TO_FASTQ -------------------------------------------------------
outs = {}
if run_log and os.path.exists(run_log):
    for line in open(run_log):
        if line.startswith("OUT\t"):
            _, sid, r1, r2, ind = line.rstrip("\n").split("\t")
            outs[sid] = (r1, r2, [i for i in ind.split(",") if i])
if not run_log:
    for case in ("illumina_dual", "illumina_single", "unpaired"):
        skip(case, "nextflow is needed")
else:
    q = os.path.join(FIX, "q", "lib_S1_L001")
    v = os.path.join(FIX, "v", "lib_S1_L001")
    if "dual" in outs:
        r1, r2, ind = outs["dual"]
        names = sorted(os.path.basename(p) for p in [r1, r2] + ind)
        want_names = sorted(f"dual_S1_L001_{r}_001.fastq.gz" for r in ("R1", "R2", "I1", "I2"))
        same = all(reads(got) == reads(f"{q}_{r}_001.fastq.gz")
                   for got, r in zip([r1, r2] + sorted(ind), ["R1", "R2", "I1", "I2"]))
        check("illumina_dual", same and names == want_names,
              f"{len(reads(r1))} pairs; {', '.join(names)}; " + ("reads identical" if same else "reads DIFFER"))
    else:
        check("illumina_dual", False, "no output, see cram_input_run.log")
    if "single" in outs:
        r1, r2, ind = outs["single"]
        ok = len(ind) == 1 and ind[0].endswith("_I1_001.fastq.gz") and reads(ind[0]) == reads(f"{v}_I2_001.fastq.gz") \
            and reads(r1) == reads(f"{v}_R1_001.fastq.gz")
        check("illumina_single", ok, f"index files: {[os.path.basename(i) for i in ind]}")
    else:
        check("illumina_single", False, "no output, see cram_input_run.log")
    import glob
    errs = "".join(open(f).read() for f in glob.glob(os.path.join(work, "nf", "work", "*", "*", ".command.err")))
    check("unpaired", "unpaired" not in outs and "not one of a read pair" in errs,
          "a single-end CRAM " + ("fails with a message" if "unpaired" not in outs else "was CONVERTED"))

# ---- Ultima CRAMs --------------------------------------------------------
d = os.path.join(FIX, "q")
lists = open(os.path.join(d, "whitelists.txt")).read().split()
subprocess.run([python, SCRIPT, "prepare", "--layout", "quantum", "--whitelists", *lists,
                "--outdir", os.path.join(d, "prep")], check=True, capture_output=True)
recs = zip(*[reads(os.path.join(d, f"lib_S1_L001_{r}_001.fastq.gz")) for r in ("R1", "R2", "I1", "I2")])
truth = {r["read"]: r for r in csv.DictReader(open(os.path.join(d, "truth.tsv")), delimiter="\t")}

def ultima(name, umi_where):
    sam = os.path.join(work, f"{name}.sam")
    rgs = set()
    body = []
    for (n, s1, q1), (_, s2, _), (_, si1, _), (_, si2, _) in zip(*[reads(os.path.join(d, f"lib_S1_L001_{r}_001.fastq.gz")) for r in ("R1", "R2", "I1", "I2")]):
        umi, rt = s2[:8], s2[8:16]
        seq, qual = (umi + s1, "I" * 8 + q1) if umi_where == "seq" else (s1, q1)
        # as the Ultima trimmer writes them ("Scale RNA" format): b2/b3 start with a constant 4 nt
        tags = [f"RG:Z:{si2}", f"b1:Z:{si1[0:8]}", f"b2:Z:CTGT{si1[12:20]}", f"b3:Z:GGTT{si1[24:32]}", f"b4:Z:{rt}"]
        if umi_where == "tag":
            tags.append(f"US:Z:{umi}")
        rgs.add(si2)
        body.append("\t".join([n, "4", "*", "0", "0", "*", "*", "0", "0", seq, qual] + tags))
    with open(sam, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n")
        fh.writelines(f"@RG\tID:{rg}\n" for rg in sorted(rgs))
        fh.write("\n".join(body) + "\n")
    cram = os.path.join(work, f"{name}.cram")
    subprocess.run(["samtools", "view", "-C", "-o", cram, sam], check=True, capture_output=True)
    out = os.path.join(work, f"{name}_out")
    p = subprocess.run([python, SCRIPT, "demux", "--barcodes", os.path.join(d, "prep", "barcodes.json"),
                        "--samples", os.path.join(d, "samples.tsv"), "--chunk", "c1", "--cram", cram,
                        "--ultima-umi", "tag:US" if umi_where == "tag" else "seq:0:8", "--outdir", out],
                       capture_output=True, text=True)
    if p.returncode != 0:
        check(name, False, p.stderr[-300:]); return
    m = json.load(open(os.path.join(out, "c1.metrics.json")))
    want = {k: v for k, v in Counter(r["category"] for r in truth.values()).items()}
    got = {k: v for k, v in m["read_categories"].items() if v}
    bad, seen = 0, 0
    for s in ("sampleA", "sampleB", "sampleC"):
        bcs = reads(os.path.join(out, f"{s}.c1_BC_UMI.fastq.gz"))
        cds = reads(os.path.join(out, f"{s}.c1_cDNA.fastq.gz"))
        for (n, bc, _), (_, cd, cq) in zip(bcs, cds):
            seen += 1
            t = truth[n]
            bad += (t["sample"], t["bc_read"], int(t["cdna_len"])) != (s, bc, len(cd)) or len(cd) != len(cq)
    check(name, got == want and seen == want["Pass"] and not bad and m["input"] == "ultima_cram",
          f"{m['reads']} records, {got}" + ("" if got == want else f" (expected {want})") + f"; {seen} written, {bad} wrong")

ultima("ultima_tag", "tag")
ultima("ultima_seq", "seq")

# ---- tag values ----------------------------------------------------------
seg = sb.Segment({"name": "pcr", "entries": [
    {"id": "QSR-1", "sequences": ["ACGTACGT", "TTGGCCAA"], "canonical": "ACGTACGT"},
    {"id": "QSR-12", "sequences": ["GGGGCCCC"], "canonical": "GGGGCCCC"}]}, 1)
cases = [("QSR-1", 0), ("QSR-12", 1), ("run7-QSR-12", 1), ("QSR-1_lane2", 0), ("ACGTACGT", 0), ("TTGGCCAT", 0),
         ("lib_GGGGCCCC", 1), ("QSR-3", None), ("NNNNNNNN", None),
         # a constant linker before (or after) the barcode, with no START given (the default tags give b2:4, b3:4)
         ("CTGTGGGGCCCC", 1), ("GGTTTTGGCCAA", 0), ("GGGGCCCATTTT", 1), ("CTGTACGTACGA", 0), ("CTGTAAAAAAAA", None),
         ("ACGTACGTGGGGCCCC", 1)]   # the end is tried first
wrong = [(v, e, seg.match_value(v)) for v, e in cases if seg.match_value(v) != e]
check("tag_values", not wrong, f"{len(cases)} tag values" + (f"; wrong: {wrong}" if wrong else " matched as expected"))
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
