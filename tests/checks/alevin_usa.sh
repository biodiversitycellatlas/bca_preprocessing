#!/usr/bin/env bash
# description: Check that alevin-fry's USA matrices are split into spliced/unspliced/ambiguous by position
#
# alevin-fry in USA mode writes three column blocks per gene: the spliced counts
# under the bare gene IDs, then <gene>-U, then <gene>-A (there is no -S). Every
# script that reads such a matrix goes through bin/alevin_usa.py; this check
# drives each of them on a fixture in exactly that layout.
#
# The failures it guards against were all silent:
#   - a collapse that did not recognise the layout copied the matrix through, so
#     every downstream step saw each gene three times;
#   - cell calling then counted every column, so --counts had no effect and a
#     gene detected as spliced and unspliced counted twice in the gene statistics;
#   - a split by suffix mistakes a gene ID ending in '-A' (HLA-A) for another
#     gene's ambiguous column. The fixture has one (gene 3).
#
# tests/lib/make_velocyto_fixture.py writes counts that identify their block,
# gene and cell, so the cases assert on values rather than shapes.
#
# No sequencing data and no alevin-fry are needed.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

COLLAPSE_TOOL="$PROJECT_ROOT/bin/collapse_alevin_usa.py"
SECONDDERIV_TOOL="$PROJECT_ROOT/bin/secondderiv_alevin.py"
FIXTURE_GEN="$TESTS_DIR/lib/make_velocyto_fixture.py"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh alevin_usa [-- OPTIONS]
       tests/checks/alevin_usa.sh [OPTIONS]

Validate the handling of alevin-fry's USA (spliced/unspliced/ambiguous) matrices.

Options:
  --python PATH   Python interpreter to use (default: \$BCA_PYTHON, else
                  python3, else python).
  --keep          Keep the generated fixtures and outputs for inspection.
  -h, --help      Show this message.

Cases:
  collapse      collapse_alevin_usa.py sums exactly the requested blocks, for
                every --counts value, and keeps the bare gene IDs
  not_usa       a matrix not in the USA layout (plain genes, or -S/-U/-A
                suffixes) fails without writing output
  cell_calling  secondderiv_alevin.py's UMI totals and filter use the requested
                blocks, keep every USA column and count each gene once
  gene_count    both dashboards report the reference's gene count, not 3x it
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

log_header "alevin_usa"

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"
    elif have_cmd python; then PYTHON="python"
    fi
fi

if [[ -z "$PYTHON" ]] || ! "$PYTHON" -c 'import sys; sys.exit(0)' >/dev/null 2>&1; then
    record SKIP "python" "no usable python interpreter found"
    finish_check
    exit 0
fi

if ! "$PYTHON" -c 'import numpy, pandas, scipy' >/dev/null 2>&1; then
    record SKIP "python" "numpy, pandas and scipy are needed to build the fixture"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/alevin_usa"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# assert NAME PYTHON_EXPR DETAIL [PATH...]
#
# PYTHON_EXPR is evaluated with the fixture's helpers in scope and must yield
# (ok, message). Paths arrive as LIB, BIN, W and then A[0], A[1], ... in the
# order given (see velocity_matrix.sh for why they are not interpolated).
assert() {
    local name="$1" expr="$2" detail="${3:-}"; shift 3 || shift $#
    local out
    if out="$("$PYTHON" -c "
import sys
LIB, BIN, W, A = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
sys.path.insert(0, LIB)
sys.path.insert(0, BIN)
import json
import numpy as np, scipy.io as sio
from make_velocyto_fixture import usa_sum, barcode, gene, N_GENES, N_CELLS

def read_lines(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]

def read_mtx(path):
    return sio.mmread(path).toarray()

def cell_total(blocks, c):
    return sum(usa_sum(blocks, g, c) for g in range(N_GENES))

ok, msg = ($expr)
print(('PASS' if ok else 'FAIL') + '\t' + str(msg))
" "$TESTS_DIR/lib" "$PROJECT_ROOT/bin" "$WORKDIR" "$@" 2>&1)"; then
        # The verdict is the last line: anything printed before it (a library warning on stderr) is not part of it
        local last="${out##*$'\n'}"
        local status="${last%%$'\t'*}" msg="${last#*$'\t'}"
        if [[ "$status" == "PASS" ]]; then
            record PASS "$name" "${detail:-$msg}"
        else
            record FAIL "$name" "$msg"
        fi
    else
        record FAIL "$name" "assertion crashed: $(printf '%s' "$out" | tail -3 | tr '\n' ' ')"
    fi
}

COUNTS_VALUES=(SUA SA S UA U)

# --------------------------------------------------------------------------
# Fixture
# --------------------------------------------------------------------------

FIXTURE="$WORKDIR/fixture"
if ! run_logged "$BCA_TEST_LOGDIR/alevin_usa_fixture.log" \
        "$PYTHON" "$FIXTURE_GEN" "$FIXTURE"; then
    record FAIL "fixture" "see alevin_usa_fixture.log"
    finish_check
    exit 1
fi
ALEVIN="$FIXTURE/alevin"

# Two matrices that are not USA sets, sharing the fixture's counts: plain gene
# columns (no blocks at all), and the -S/-U/-A naming no alevin-fry writes
"$PYTHON" - "$ALEVIN" "$WORKDIR" <<'PY'
import os, shutil, sys
src, work = sys.argv[1], sys.argv[2]
cols = [l.strip() for l in open(os.path.join(src, "quants_mat_cols.txt")) if l.strip()]
n = len(cols) // 3
for name, names in (("plain", [f"P{i}" for i in range(len(cols))]),
                    ("suffix_s", [g + "-S" for g in cols[:n]] + cols[n:])):
    out = os.path.join(work, name)
    os.makedirs(out, exist_ok=True)
    for f in ("quants_mat.mtx", "quants_mat_rows.txt"):
        shutil.copyfile(os.path.join(src, f), os.path.join(out, f))
    with open(os.path.join(out, "quants_mat_cols.txt"), "w") as fh:
        fh.write("\n".join(names) + "\n")
PY
record PASS "fixture" "USA matrix (bare / -U / -A) and two non-USA layouts written"

# --------------------------------------------------------------------------
# Case: collapse
# --------------------------------------------------------------------------

for counts in "${COUNTS_VALUES[@]}"; do
    out="$WORKDIR/collapse_$counts"
    if run_logged "$BCA_TEST_LOGDIR/alevin_usa_collapse_$counts.log" \
            "$PYTHON" "$COLLAPSE_TOOL" --dir "$ALEVIN" --outdir "$out" --counts "$counts"; then
        assert "collapse.$counts" \
            "(read_lines(A[0] + '/quants_mat_cols.txt') == [gene(i) for i in range(N_GENES)]
              and read_lines(A[0] + '/quants_mat_rows.txt') == [barcode(c) for c in range(N_CELLS)]
              and all(read_mtx(A[0] + '/quants_mat.mtx')[c][g] == usa_sum(A[1], g, c)
                      for g in range(N_GENES) for c in range(N_CELLS)),
              f\"cols={read_lines(A[0] + '/quants_mat_cols.txt')}\")" \
            "one column per gene holding exactly $counts, under the bare gene IDs" \
            "$out" "$counts"
    else
        record FAIL "collapse.$counts" "see alevin_usa_collapse_$counts.log"
    fi
done

# --------------------------------------------------------------------------
# Case: not_usa
#
# The collapse used to copy such a matrix through with a warning and exit 0.
# --------------------------------------------------------------------------

for layout in plain suffix_s; do
    out="$WORKDIR/not_usa_$layout"
    if run_logged "$BCA_TEST_LOGDIR/alevin_usa_not_usa_$layout.log" \
            "$PYTHON" "$COLLAPSE_TOOL" --dir "$WORKDIR/$layout" --outdir "$out" --counts SUA; then
        record FAIL "not_usa.collapse_$layout" "a non-USA matrix was accepted"
    elif [[ -e "$out/quants_mat.mtx" ]]; then
        record FAIL "not_usa.collapse_$layout" "failed but still wrote a matrix"
    else
        record PASS "not_usa.collapse_$layout" "rejected without output"
    fi

    if run_logged "$BCA_TEST_LOGDIR/alevin_usa_not_usa_umis_$layout.log" \
            "$PYTHON" "$SECONDDERIV_TOOL" umis --dir "$WORKDIR/$layout" \
            --output "$WORKDIR/not_usa_umis_$layout.txt" --counts SUA; then
        record FAIL "not_usa.cell_calling_$layout" "a non-USA matrix was accepted"
    else
        record PASS "not_usa.cell_calling_$layout" "rejected"
    fi
done

# --------------------------------------------------------------------------
# Case: cell_calling
#
# Totals per cell rise with the cell index in the fixture, so a cutoff at the
# total of cell 5 must keep exactly cells 5, 6 and 7, for every --counts.
# --------------------------------------------------------------------------

for counts in "${COUNTS_VALUES[@]}"; do
    umis="$WORKDIR/umis_$counts.txt"
    if run_logged "$BCA_TEST_LOGDIR/alevin_usa_umis_$counts.log" \
            "$PYTHON" "$SECONDDERIV_TOOL" umis --dir "$ALEVIN" --output "$umis" --counts "$counts"; then
        assert "cell_calling.umis_$counts" \
            "([int(v) for v in read_lines(A[0])] == sorted((cell_total(A[1], c) for c in range(N_CELLS)), reverse=True),
              read_lines(A[0]))" \
            "per-cell $counts totals, descending" \
            "$umis" "$counts"
    else
        record FAIL "cell_calling.umis_$counts" "see alevin_usa_umis_$counts.log"
        continue
    fi

    cutoff="$("$PYTHON" -c "
import sys; sys.path.insert(0, sys.argv[1])
from make_velocyto_fixture import usa_sum, N_GENES
print(sum(usa_sum(sys.argv[2], g, 5) for g in range(N_GENES)))
" "$TESTS_DIR/lib" "$counts")"

    out="$WORKDIR/filter_$counts"
    stats="$WORKDIR/filter_$counts.json"
    if run_logged "$BCA_TEST_LOGDIR/alevin_usa_filter_$counts.log" \
            "$PYTHON" "$SECONDDERIV_TOOL" filter --dir "$ALEVIN" --cutoff "$cutoff" \
            --outdir "$out" --stats "$stats" --counts "$counts"; then
        assert "cell_calling.filter_$counts" \
            "(read_lines(A[0] + '/quants_mat_rows.txt') == [barcode(c) for c in (5, 6, 7)]
              and read_lines(A[0] + '/quants_mat_cols.txt') == read_lines(A[2] + '/quants_mat_cols.txt')
              and read_mtx(A[0] + '/quants_mat.mtx').shape == (3, 3 * N_GENES),
              read_lines(A[0] + '/quants_mat_rows.txt'))" \
            "cells 5-7 kept, every USA column kept" \
            "$out" "$counts" "$ALEVIN"

        assert "cell_calling.stats_$counts" \
            "(json.load(open(A[0]))['total_genes_detected'] == N_GENES
              and json.load(open(A[0]))['median_genes_per_cell'] == N_GENES
              and json.load(open(A[0]))['estimated_cells'] == 3
              and json.load(open(A[0]))['median_umis_per_cell'] == cell_total(A[1], 6),
              json.load(open(A[0])))" \
            "genes counted once, UMIs on the $counts basis" \
            "$stats" "$counts"
    else
        record FAIL "cell_calling.filter_$counts" "see alevin_usa_filter_$counts.log"
    fi
done

# --------------------------------------------------------------------------
# Case: gene_count
# --------------------------------------------------------------------------

assert "gene_count" \
    "(lambda gd, ms, P: (
        (gd.count_alevin_genes({}, A[0]), ms.count_alevin_genes({}, P(A[0])),
         gd.count_alevin_genes({}, A[1]), ms.count_alevin_genes({}, P(A[1])),
         gd.count_alevin_genes({'num_genes': 3 * N_GENES, 'usa_mode': True}, None),
         ms.count_alevin_genes({'num_genes': 3 * N_GENES, 'usa_mode': True}, None))
        == (N_GENES, N_GENES, 3 * N_GENES, 3 * N_GENES, N_GENES, N_GENES),
        'USA columns, plain columns, quant.json fallback'))(
        __import__('generate_dashboard'), __import__('dashboard_mappingstats'), __import__('pathlib').Path)" \
    "both dashboards: N genes for a USA set, every name otherwise, num_genes/3 from quant.json" \
    "$ALEVIN/quants_mat_cols.txt" "$WORKDIR/plain/quants_mat_cols.txt"

finish_check
