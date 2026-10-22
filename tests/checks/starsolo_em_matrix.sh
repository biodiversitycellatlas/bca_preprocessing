#!/usr/bin/env bash
# description: Check that STARsolo's multi-gene matrix (UniqueAndMult-EM.mtx) is the one curve, cutoff and filter use
#
# With star_soloMultiMappers = "EM" (every protocol's default) STAR writes two raw matrices:
# matrix.mtx with the unique-gene counts and UniqueAndMult-EM.mtx with the unique plus the
# multi-gene counts, fractional. The pipeline analyses the second (lib/StarsoloMatrix.groovy),
# so the UMI curve the cutoff is read off, the filter and the subset to STAR's own cells
# must all read it, and on the same basis. The failures guarded against are silent:
#   - a curve from STAR's UMIperCellSorted.txt (unique counts) with a filter on EM totals
#     keeps a different cell set than the curve shows;
#   - "genes detected" counting every EM fraction above zero inflates the gene counts;
#   - a missing multi-gene matrix silently falling back to the unique one;
#   - with STAR's own cell call, the dashboard reporting the median UMIs and genes of the
#     unique counts (Summary.csv) for the EM matrix that is analysed.
#
# The fixture has a cell whose total clears the cutoff only with its EM counts, and a gene
# that only ever receives EM fractions below one UMI.
#
# Needs numpy, pandas and scipy. The name case also needs nextflow, and is skipped without it.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

FILTER_TOOL="$PROJECT_ROOT/bin/secondderiv_filter_matrices.py"
SUBSET_TOOL="$PROJECT_ROOT/bin/subset_matrices_to_cells.py"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh starsolo_em_matrix [-- OPTIONS]
       tests/checks/starsolo_em_matrix.sh [OPTIONS]

Validate the use of STARsolo's multi-gene matrix in cell calling.

Options:
  --python PATH   Python interpreter to use (default: \$BCA_PYTHON, else
                  python3, else python).
  --keep          Keep the generated fixture and outputs for inspection.
  -h, --help      Show this message.

Cases:
  umis            secondderiv_filter_matrices.py umis writes the descending totals
                  of the requested matrix, truncated to whole UMIs
  filter          the cutoff keeps the cells whose EM total reaches it, the
                  filtered matrix holds the EM counts as matrix.mtx, and genes
                  count as detected from one UMI
  missing_matrix  asking for a multi-gene matrix that is not there fails
  subset          STAR's own cells, subset from the EM matrix as
                  SUBSET_STARSOLO_CELLS does, and that matrix's statistics
  dashboard       dashboard_mappingstats.py takes the median UMIs, median genes
                  and genes detected from those statistics, the rest from STAR
  name           StarsoloMatrix picks the file from star_soloMultiMappers
                  (needs nextflow)
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

log_header "starsolo_em_matrix"

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"
    elif have_cmd python; then PYTHON="python"
    fi
fi

if [[ -z "$PYTHON" ]] || ! "$PYTHON" -c 'import numpy, pandas, scipy' >/dev/null 2>&1; then
    record SKIP "python" "a python with numpy, pandas and scipy is needed"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/starsolo_em_matrix"
rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# assert NAME DETAIL PYTHON_CODE: the code runs with W (the work directory) and the fixture
# in scope, and sets ok and msg
assert() {
    local name="$1" detail="$2" code="$3" out
    if out="$("$PYTHON" - "$WORKDIR" 2>&1 <<PY
import json, os, sys
import numpy as np, scipy.io as sio
W = sys.argv[1]
fx = json.load(open(os.path.join(W, "fixture.json")))
U, E = np.array(fx["unique"]), np.array(fx["em"])           # genes x cells
BC, CUT = fx["barcodes"], fx["cutoff"]
def lines(p): return [l.strip() for l in open(p) if l.strip()]
def mtx(p): return sio.mmread(p).toarray()
ok, msg = False, "no verdict"
$code
print(("PASS" if ok else "FAIL") + "\t" + str(msg))
PY
)"; then
        local last="${out##*$'\n'}"
        local status="${last%%$'\t'*}" msg="${last#*$'\t'}"
        if [[ "$status" == "PASS" ]]; then record PASS "$name" "${detail:-$msg}"; else record FAIL "$name" "$msg"; fi
    else
        record FAIL "$name" "assertion crashed: $(printf '%s' "$out" | tail -3 | tr '\n' ' ')"
    fi
}

# --------------------------------------------------------------------------
# Fixture: a STARsolo raw/ directory with both matrices
# --------------------------------------------------------------------------

RAW="$WORKDIR/raw"
"$PYTHON" - "$WORKDIR" <<'PY'
import json, os, sys
import numpy as np, scipy.io as sio, scipy.sparse as sp
W = sys.argv[1]
raw = os.path.join(W, "raw"); os.makedirs(raw)
barcodes = [f"CELL{i:02d}" for i in range(8)]
genes = [f"G{i}" for i in range(5)]
# genes x cells. Cell 5 has 9 unique UMIs and reaches the cutoff of 10 only with EM;
# gene G4 only ever receives fractions below one UMI.
unique = np.array([
    [9, 7, 5, 4, 1, 3, 0, 0],
    [6, 5, 4, 3, 1, 3, 1, 0],
    [4, 3, 3, 2, 0, 3, 0, 1],
    [1, 0, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
], dtype=float)
extra = np.zeros_like(unique)
extra[0, 5], extra[3, 5] = 0.75, 0.75      # cell 5: 9 -> 10.5
extra[3, 0], extra[3, 2] = 2.5, 1.5
extra[4, 0], extra[4, 1], extra[4, 5] = 0.4, 0.3, 0.2
em = unique + extra
sio.mmwrite(os.path.join(raw, "matrix.mtx"), sp.coo_matrix(unique.astype(int)))
sio.mmwrite(os.path.join(raw, "UniqueAndMult-EM.mtx"), sp.coo_matrix(em))
with open(os.path.join(raw, "barcodes.tsv"), "w") as fh: fh.write("\n".join(barcodes) + "\n")
with open(os.path.join(raw, "features.tsv"), "w") as fh:
    fh.write("".join(f"{g}\t{g}\tGene Expression\n" for g in genes))
# STAR's own call, as filtered/barcodes.tsv would list it (out of matrix order on purpose)
with open(os.path.join(W, "star_cells.tsv"), "w") as fh: fh.write("CELL03\nCELL00\nCELL05\n")
# The unique-only directory, for the missing-matrix case
os.makedirs(os.path.join(W, "raw_unique"))
for f in ("matrix.mtx", "barcodes.tsv", "features.tsv"):
    os.link(os.path.join(raw, f), os.path.join(W, "raw_unique", f))
json.dump({"unique": unique.tolist(), "em": em.tolist(), "barcodes": barcodes, "cutoff": 10},
          open(os.path.join(W, "fixture.json"), "w"))
PY
record PASS "fixture" "raw/ with matrix.mtx and UniqueAndMult-EM.mtx (8 barcodes x 5 genes)"

# --------------------------------------------------------------------------
# Case: umis
# --------------------------------------------------------------------------

if run_logged "$BCA_TEST_LOGDIR/starsolo_em_umis.log" \
        "$PYTHON" "$FILTER_TOOL" umis -d "$RAW" -m UniqueAndMult-EM.mtx -o "$WORKDIR/em_curve.txt" \
   && run_logged "$BCA_TEST_LOGDIR/starsolo_em_umis_unique.log" \
        "$PYTHON" "$FILTER_TOOL" umis -d "$RAW" -o "$WORKDIR/unique_curve.txt"; then
    assert "umis" "EM curve = sorted floor(EM totals); default matrix.mtx = unique totals" '
em = [int(x) for x in lines(os.path.join(W, "em_curve.txt"))]
un = [int(x) for x in lines(os.path.join(W, "unique_curve.txt"))]
want_em = sorted(np.floor(E.sum(axis=0)).astype(int).tolist(), reverse=True)
want_un = sorted(U.sum(axis=0).astype(int).tolist(), reverse=True)
ok = em == want_em and un == want_un
msg = f"em={em} want {want_em}; unique={un} want {want_un}"
'
else
    record FAIL "umis" "see starsolo_em_umis*.log"
fi

# --------------------------------------------------------------------------
# Case: filter
# --------------------------------------------------------------------------

if run_logged "$BCA_TEST_LOGDIR/starsolo_em_filter.log" \
        "$PYTHON" "$FILTER_TOOL" filter -d "$RAW" -m UniqueAndMult-EM.mtx -c 10 \
            -o "$WORKDIR/filtered_em" -s "$WORKDIR/em_stats.json" \
   && run_logged "$BCA_TEST_LOGDIR/starsolo_em_filter_unique.log" \
        "$PYTHON" "$FILTER_TOOL" filter -d "$RAW" -c 10 \
            -o "$WORKDIR/filtered_unique" -s "$WORKDIR/unique_stats.json"; then
    assert "filter.cells" "cell 5 kept on its EM total (10.5), dropped on its unique total (9)" '
keep_em = [BC[i] for i in np.where(E.sum(axis=0) >= CUT)[0]]
keep_un = [BC[i] for i in np.where(U.sum(axis=0) >= CUT)[0]]
got_em = lines(os.path.join(W, "filtered_em", "barcodes.tsv"))
got_un = lines(os.path.join(W, "filtered_unique", "barcodes.tsv"))
ok = got_em == keep_em and got_un == keep_un and "CELL05" in got_em and "CELL05" not in got_un
msg = f"em kept {got_em}, unique kept {got_un}"
'
    assert "filter.matrix" "written as matrix.mtx, real-valued, holding the EM counts" '
p = os.path.join(W, "filtered_em", "matrix.mtx")
keep = np.where(E.sum(axis=0) >= CUT)[0]
header = open(p).readline()
ok = "real" in header and np.allclose(mtx(p), E[:, keep]) and not os.path.exists(os.path.join(W, "filtered_em", "UniqueAndMult-EM.mtx"))
msg = f"header {header.strip()!r}"
'
    assert "filter.genes_detected" "genes count from one UMI: G4 (fractions only) is never detected" '
st = json.load(open(os.path.join(W, "em_stats.json")))
keep = np.where(E.sum(axis=0) >= CUT)[0]
want_median = float(np.median((E[:, keep] >= 1).sum(axis=0)))
want_total = int((E[:, keep].sum(axis=1) >= 1).sum())
ok = st["median_genes_per_cell"] == want_median and st["total_genes_detected"] == want_total == 4 \
     and st["median_umis_per_cell"] == float(np.median(E[:, keep].sum(axis=0))) \
     and st["umi_threshold_applied"] == CUT and st["estimated_cells"] == len(keep)
msg = f"stats {st}"
'
else
    record FAIL "filter" "see starsolo_em_filter*.log"
fi

# --------------------------------------------------------------------------
# Case: missing_matrix
# --------------------------------------------------------------------------

if "$PYTHON" "$FILTER_TOOL" filter -d "$WORKDIR/raw_unique" -m UniqueAndMult-EM.mtx -c 10 \
        -o "$WORKDIR/filtered_missing" -s "$WORKDIR/missing_stats.json" >"$BCA_TEST_LOGDIR/starsolo_em_missing.log" 2>&1; then
    record FAIL "missing_matrix" "a missing UniqueAndMult-EM.mtx did not fail"
elif [[ -e "$WORKDIR/filtered_missing/matrix.mtx" ]]; then
    record FAIL "missing_matrix" "failed, but still wrote a matrix"
elif grep -q "star_soloMultiMappers" "$BCA_TEST_LOGDIR/starsolo_em_missing.log"; then
    record PASS "missing_matrix" "fails, naming star_soloMultiMappers, without writing output"
else
    record FAIL "missing_matrix" "failed without naming the setting: $(tail -1 "$BCA_TEST_LOGDIR/starsolo_em_missing.log")"
fi

# --------------------------------------------------------------------------
# Case: subset (the commands SUBSET_STARSOLO_CELLS runs)
# --------------------------------------------------------------------------

if run_logged "$BCA_TEST_LOGDIR/starsolo_em_subset.log" \
        "$PYTHON" "$SUBSET_TOOL" --dir "$RAW" --barcodes "$WORKDIR/star_cells.tsv" \
            --matrices UniqueAndMult-EM.mtx --outdir "$WORKDIR/subset" \
   && mv "$WORKDIR/subset/UniqueAndMult-EM.mtx" "$WORKDIR/subset/matrix.mtx" \
   && run_logged "$BCA_TEST_LOGDIR/starsolo_em_subset_stats.log" \
        "$PYTHON" "$FILTER_TOOL" stats -d "$WORKDIR/subset" -m matrix.mtx --label UniqueAndMult-EM.mtx \
            -s "$WORKDIR/S1_UniqueAndMult-EM_statistics.json"; then
    assert "subset" "STAR's cells, in its order, from the EM matrix only" '
order = [BC.index(b) for b in ["CELL03", "CELL00", "CELL05"]]
got = lines(os.path.join(W, "subset", "barcodes.tsv"))
files = sorted(f for f in os.listdir(os.path.join(W, "subset")) if f.endswith(".mtx"))
ok = got == ["CELL03", "CELL00", "CELL05"] and files == ["matrix.mtx"] \
     and np.allclose(mtx(os.path.join(W, "subset", "matrix.mtx")), E[:, order])
msg = f"barcodes {got}, files {files}"
'
    assert "subset.stats" "EM statistics of STAR's cells, labelled with the raw matrix, no cutoff" '
st = json.load(open(os.path.join(W, "S1_UniqueAndMult-EM_statistics.json")))
S = E[:, [BC.index(b) for b in ["CELL03", "CELL00", "CELL05"]]]
want = {"matrix": "UniqueAndMult-EM.mtx", "estimated_cells": 3,
        "mean_umis_per_cell": float(np.mean(S.sum(axis=0))),
        "median_umis_per_cell": float(np.median(S.sum(axis=0))),
        "median_genes_per_cell": float(np.median((S >= 1).sum(axis=0))),
        "total_genes_detected": int((S.sum(axis=1) >= 1).sum())}
ok = st.keys() == want.keys() and all(np.isclose(st[k], v) if k != "matrix" else st[k] == v for k, v in want.items())
msg = f"stats {st}, want {want}"
'
else
    record FAIL "subset" "see starsolo_em_subset*.log"
fi

# --------------------------------------------------------------------------
# Case: dashboard (dashboard_mappingstats.py under star_solocellfilter)
# --------------------------------------------------------------------------

# An outdir as MAPPING_STATS stages it: Summary.csv with STAR's unique-count figures,
# filtered/ (STAR's cells, here the subset), Log.out with nUMImin and the statistics JSON
DASH="$WORKDIR/dash"
GENE="$DASH/mapping_STARsolo/S1/S1_Solo.out/GeneFull_Ex50pAS"
if [[ -s "$WORKDIR/S1_UniqueAndMult-EM_statistics.json" ]]; then
    mkdir -p "$GENE"
    cp -r "$WORKDIR/subset" "$GENE/filtered"
    cp "$WORKDIR/S1_UniqueAndMult-EM_statistics.json" "$GENE/"
    cat >"$GENE/Summary.csv" <<'CSV'
Number of Reads,1000
Sequencing Saturation,0.25
Estimated Number of Cells,3
Mean Reads per Cell,300
Median UMI per Cell,1
Median GeneFull_Ex50pAS per Cell,1
Total GeneFull_Ex50pAS Detected,1
CSV
    printf 'Starting Solo post-map for GeneFull_Ex50pAS\ncellFiltering: nUMImin=7\n' >"$DASH/mapping_STARsolo/S1/S1_Log.out"
fi

if [[ ! -d "$GENE" ]]; then
    record SKIP "dashboard" "needs the subset case's statistics"
elif (cd "$WORKDIR" && run_logged "$BCA_TEST_LOGDIR/starsolo_em_dashboard.log" \
        "$PYTHON" "$PROJECT_ROOT/bin/dashboard_mappingstats.py" "$DASH" -o "$WORKDIR/mapping_stats.tsv"); then
    assert "dashboard" "median UMIs, median genes and genes detected from the EM statistics; cells and saturation from STAR" '
import pandas as pd
row = pd.read_csv(os.path.join(W, "mapping_stats.tsv"), sep="\t", dtype=str).set_index("Sample").loc["S1"]
st = json.load(open(os.path.join(W, "dash/mapping_STARsolo/S1/S1_Solo.out/GeneFull_Ex50pAS/S1_UniqueAndMult-EM_statistics.json")))
got = {c: row[c] for c in ("Median UMI Counts per Cell", "Median Genes per Cell", "Total Genes Detected", "N cells", "Saturation", "Mean Reads per Cell")}
want = {"Median UMI Counts per Cell": str(int(round(st["median_umis_per_cell"]))),
        "Median Genes per Cell": str(int(round(st["median_genes_per_cell"]))),
        "Total Genes Detected": str(st["total_genes_detected"]),
        "N cells": "3", "Mean Reads per Cell": "300"}
ok = all(got[k] == v for k, v in want.items()) and got["Saturation"].startswith("25")
msg = f"got {got}, want {want} and Saturation 25%"
'
else
    record FAIL "dashboard" "see starsolo_em_dashboard.log"
fi

# --------------------------------------------------------------------------
# Case: name (lib/StarsoloMatrix.groovy, through a minimal Nextflow project)
# --------------------------------------------------------------------------

if ! have_cmd nextflow; then
    record SKIP "name" "nextflow not found"
else
    NF="$WORKDIR/nf"
    mkdir -p "$NF"
    ln -s "$PROJECT_ROOT/lib" "$NF/lib"
    cat >"$NF/main.nf" <<'NF'
def proto = [p10: [star_soloMultiMappers: 'EM'], pU: [star_soloMultiMappers: 'Unique']]
def cases = [
    [[star_soloMultiMappers: null,          protocol: 'p10', cellfilter_method: 'second_derivative'],   'UniqueAndMult-EM.mtx', true],
    [[star_soloMultiMappers: 'Unique',      protocol: 'p10', cellfilter_method: 'second_derivative'],   'matrix.mtx',           false],
    [[star_soloMultiMappers: null,          protocol: 'pU',  cellfilter_method: 'second_derivative'],   'matrix.mtx',           false],
    [[star_soloMultiMappers: 'Unique Uniform EM', protocol: 'pU', cellfilter_method: 'manual_cutoff'], 'UniqueAndMult-Uniform.mtx', true],
    [[star_soloMultiMappers: 'EM',          protocol: 'pU',  cellfilter_method: 'star_solocellfilter'], 'UniqueAndMult-EM.mtx', false],
    [[star_soloMultiMappers: null,          protocol: 'none', cellfilter_method: 'second_derivative'],  'matrix.mtx',           false],
]
cases.each { p, file, curve ->
    def params_map = p + [seqtech_parameters: proto]
    def got_file = StarsoloMatrix.rawFileName(params_map)
    def got_curve = StarsoloMatrix.curveFromMatrix(params_map)
    println "${got_file == file && got_curve == curve ? 'OK' : 'BAD'}\t${p}\t${got_file}\t${got_curve}"
}
NF
    if (cd "$NF" && run_logged "$BCA_TEST_LOGDIR/starsolo_em_name.log" nextflow -q run main.nf); then
        if grep -q "^BAD" "$BCA_TEST_LOGDIR/starsolo_em_name.log"; then
            record FAIL "name" "$(grep '^BAD' "$BCA_TEST_LOGDIR/starsolo_em_name.log" | head -2 | tr '\n' ' ')"
        elif [[ "$(grep -c '^OK' "$BCA_TEST_LOGDIR/starsolo_em_name.log")" -eq 6 ]]; then
            record PASS "name" "file and curve source follow star_soloMultiMappers and cellfilter_method (6 settings)"
        else
            record FAIL "name" "expected 6 results, see starsolo_em_name.log"
        fi
    else
        record FAIL "name" "nextflow run failed, see starsolo_em_name.log"
    fi
fi

finish_check
