#!/usr/bin/env bash
# description: Check that copies of reference sequence are dropped from ref_fasta_addfeature and the rRNA regions are right
#
# Rfam-derived rRNA references annotate loci of the genome itself (gene rows, no
# transcript_id) and carry the same loci, cut out of the same assembly, as the
# added FASTA. Appended as they are, every read on those loci multimaps.
# prepare_addfeature.py (PREPARE_ADDFEATURE) drops such copies and writes the
# rRNA regions CALC_READ_METRICS counts. The fixture holds one record of every
# kind, so each decision is asserted by name:
#
#   SSU_plus    exact copy of chrA:101-400 (+), named by the GTF     dropped
#   LSU_minus   exact copy of chrA:1001-1300 (-), reverse complement dropped
#   5S_mut      chrA:1501-1600 with one base changed                 kept, warned
#   rDNA        sequence absent from the genome                      kept
#   rDNA2       identical to rDNA                                    dropped
#   chrM        identical to the base contig of that name            dropped, no region
#   chrM_copy   identical to chrM under another name                 dropped
#
# chrA:301-400 is soft-masked (lower case), so the comparison must ignore case.
# The GTF also has a row on chrM with no record (a locus only), one on the added
# contig rDNA, and one on a contig of neither FASTA (ignored, warned).

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

TOOL="$PROJECT_ROOT/bin/prepare_addfeature.py"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh addfeature_prepare [-- OPTIONS]
       tests/checks/addfeature_prepare.sh [OPTIONS]

Validate PREPARE_ADDFEATURE's duplicate detection and rRNA regions on a
hand-written genome and Rfam-style added reference.

Options:
  --python PATH   Python interpreter (default: \$BCA_PYTHON, else python3).
  --keep          Keep the generated fixture and outputs for inspection.
  -h, --help      Show this message.

Cases:
  fasta           only records with no copy in the reference are kept
  report          every record's decision and the ignored GTF contig
  regions         loci, added contigs and dropped copies, 0-based
  warnings        the mutated copy and the unknown contig are warned about
  gtf_only        without an added FASTA, only the loci are written
  name_clash      a record named like a base contig but differing stops the run
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

log_header "addfeature_prepare"

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"; elif have_cmd python; then PYTHON="python"; fi
fi
if [[ -z "$PYTHON" ]]; then
    record SKIP "tools" "a Python is needed"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/addfeature_prepare"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# --------------------------------------------------------------------------
# Fixture
# --------------------------------------------------------------------------

FIX="$WORKDIR/fixture"
mkdir -p "$FIX"

"$PYTHON" - "$FIX" <<'PY'
import random, sys
fix = sys.argv[1]
rng = random.Random(1)
seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

chrA, chrM, rdna = seq(2000), seq(600), seq(500)
masked = chrA[:300] + chrA[300:400].lower() + chrA[400:]

def fasta(path, records):
    with open(path, "w") as fh:
        for name, s in records:
            fh.write(f">{name} description\n")
            fh.write("".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))

mutated = chrA[1500:1600]
mutated = mutated[:50] + ("A" if mutated[50] != "A" else "C") + mutated[51:]
fasta(f"{fix}/base.fa", [("chrA", masked), ("chrM", chrM)])
fasta(f"{fix}/add.fa", [
    ("SSU_plus", chrA[100:400]),
    ("LSU_minus", rc(chrA[1000:1300])),
    ("5S_mut", mutated),
    ("rDNA", rdna),
    ("rDNA2", rdna),
    ("chrM", chrM),
    ("chrM_copy", chrM),
])
fasta(f"{fix}/clash.fa", [("chrA", seq(300))])

rows = [
    ("chrA", 101, 400, "+", "SSU_plus"),
    ("chrA", 1001, 1300, "-", "LSU_minus"),
    ("chrA", 1501, 1600, "+", "5S_mut"),
    ("chrM", 51, 300, "+", "mt_LSU"),
    ("rDNA", 1, 500, "+", "rDNA_ann"),
    ("scaffold_x", 1, 100, "+", "5S_unknown"),
]
with open(f"{fix}/add.gtf", "w") as fh:
    for contig, start, end, strand, gid in rows:
        fh.write(f'{contig}\tRfam\tgene\t{start}\t{end}\t.\t{strand}\t.\t'
                 f'gene_id "{gid}"; gene_biotype "rRNA";\n')
PY

# run_prepare NAME OUTDIR ARGS... -- runs the tool in its own directory
run_prepare() {
    local name="$1" out="$2"; shift 2
    mkdir -p "$out"
    (cd "$out" && run_logged "$BCA_TEST_LOGDIR/addfeature_prepare_${name}.log" \
        "$PYTHON" "$TOOL" --base-fasta "$FIX/base.fa" "$@" \
            --out-fasta addfeature.fasta --out-regions rrna_regions.tsv \
            --out-report addfeature_report.tsv)
}

# check NAME PYTHON-SNIPPET OUTDIR -- PASS when the snippet prints "ok"
check() {
    local name="$1" snippet="$2" out="$3" result
    if result="$("$PYTHON" - "$out" <<<"$snippet")" && [[ "$result" == "ok" ]]; then
        record PASS "$name"
    else
        record FAIL "$name" "$result"
    fi
}

# --------------------------------------------------------------------------
# Cases: fasta, report, regions, warnings
# --------------------------------------------------------------------------

OUT="$WORKDIR/both"
if run_prepare both "$OUT" --add-fasta "$FIX/add.fa" --add-gtf "$FIX/add.gtf"; then
    check fasta '
import sys
names = [l[1:].strip() for l in open(sys.argv[1] + "/addfeature.fasta") if l.startswith(">")]
print("ok" if names == ["5S_mut", "rDNA"] else f"kept {names}, want [5S_mut, rDNA]")
' "$OUT"

    check report '
import csv, sys
rows = {(r["item"], r["name"]): r for r in csv.DictReader(open(sys.argv[1] + "/addfeature_report.tsv"), delimiter="\t")}
want = {
    ("fasta_record", "SSU_plus"): ("duplicate", "chrA:101-400(+)"),
    ("fasta_record", "LSU_minus"): ("duplicate", "chrA:1001-1300(-)"),
    ("fasta_record", "5S_mut"): ("kept", "differs from chrA:1501-1600"),
    ("fasta_record", "rDNA"): ("kept", "no copy"),
    ("fasta_record", "rDNA2"): ("duplicate", "added record rDNA"),
    ("fasta_record", "chrM"): ("duplicate", "ref_fasta contig chrM"),
    ("fasta_record", "chrM_copy"): ("duplicate", "ref_fasta contig chrM"),
    ("gtf_contig", "scaffold_x"): ("dropped", "1 row"),
}
bad = []
for key, (status, detail) in want.items():
    got = rows.get(key, {"status": None, "detail": ""})
    if got["status"] != status or detail not in got["detail"]:
        bad.append(key[1] + ": " + str(got["status"]) + " / " + got["detail"])
if set(rows) != set(want):
    bad.append(f"rows {sorted(rows)}")
print("; ".join(bad) if bad else "ok")
' "$OUT"

    # rDNA_ann stays as a locus of a kept added contig; chrM is the base contig itself,
    # so its copy adds no region, and scaffold_x is ignored
    check regions '
import csv, sys
rows = sorted((r["contig"], int(r["start"]), int(r["end"]), r["name"], r["kind"])
              for r in csv.DictReader(open(sys.argv[1] + "/rrna_regions.tsv"), delimiter="\t"))
want = sorted([
    ("chrA", 100, 400, "SSU_plus", "feature"),
    ("chrA", 1000, 1300, "LSU_minus", "feature"),
    ("chrA", 1500, 1600, "5S_mut", "feature"),
    ("chrM", 50, 300, "mt_LSU", "feature"),
    ("rDNA", 0, 500, "rDNA_ann", "feature"),
    ("5S_mut", 0, 100, "5S_mut", "added_contig"),
    ("rDNA", 0, 500, "rDNA", "added_contig"),
    ("SSU_plus", 0, 300, "SSU_plus", "duplicate_contig"),
    ("LSU_minus", 0, 300, "LSU_minus", "duplicate_contig"),
    ("rDNA2", 0, 500, "rDNA2", "duplicate_contig"),
    ("chrM_copy", 0, 600, "chrM_copy", "duplicate_contig"),
])
extra = [r for r in rows if r not in want]
missing = [r for r in want if r not in rows]
print("ok" if not extra and not missing else f"extra {extra}; missing {missing}")
' "$OUT"

    LOG="$BCA_TEST_LOGDIR/addfeature_prepare_both.log"
    if grep -q "added record 5S_mut differs" "$LOG" \
            && grep -q "contigs of neither FASTA and are ignored: scaffold_x" "$LOG" \
            && grep -q "5 record(s) of ref_fasta_addfeature copy reference sequence" "$LOG"; then
        record PASS "warnings" "mutated copy, unknown contig and dropped copies reported"
    else
        record FAIL "warnings" "expected warnings missing, see addfeature_prepare_both.log"
        tail_log "$LOG"
    fi
else
    record FAIL "both" "prepare_addfeature.py failed, see addfeature_prepare_both.log"
    tail_log "$BCA_TEST_LOGDIR/addfeature_prepare_both.log"
fi

# --------------------------------------------------------------------------
# Cases: gtf_only, name_clash
# --------------------------------------------------------------------------

OUT="$WORKDIR/gtf_only"
if run_prepare gtf_only "$OUT" --add-gtf "$FIX/add.gtf"; then
    # Without the FASTA, rDNA is a contig of neither FASTA and its row is ignored too
    check gtf_only '
import csv, os, sys
rows = sorted((r["contig"], r["name"], r["kind"])
              for r in csv.DictReader(open(sys.argv[1] + "/rrna_regions.tsv"), delimiter="\t"))
want = sorted([("chrA", "SSU_plus", "feature"), ("chrA", "LSU_minus", "feature"),
               ("chrA", "5S_mut", "feature"), ("chrM", "mt_LSU", "feature")])
empty = os.path.getsize(sys.argv[1] + "/addfeature.fasta") == 0
print("ok" if rows == want and empty else f"regions {rows}; fasta empty: {empty}")
' "$OUT"
else
    record FAIL "gtf_only" "prepare_addfeature.py failed, see addfeature_prepare_gtf_only.log"
    tail_log "$BCA_TEST_LOGDIR/addfeature_prepare_gtf_only.log"
fi

if run_prepare name_clash "$WORKDIR/name_clash" --add-fasta "$FIX/clash.fa"; then
    record FAIL "name_clash" "a record named chrA with another sequence was accepted"
elif grep -q "added record 'chrA' has the name of a ref_fasta contig" \
        "$BCA_TEST_LOGDIR/addfeature_prepare_name_clash.log"; then
    record PASS "name_clash" "stopped with a message naming the record"
else
    record FAIL "name_clash" "failed without the expected message"
    tail_log "$BCA_TEST_LOGDIR/addfeature_prepare_name_clash.log"
fi

finish_check
