#!/usr/bin/env bash
# description: Check that the filtering chain's triplets and the annotated h5ad it ends in are packed correctly
#
# MTX_TO_H5AD now runs last: doublet detection, doublet filtering and CellSweep
# all read and write the (mtx, barcodes, features) triplet, and their results are
# merged into one AnnData object at the end. Every failure mode along that chain
# is silent. STARsolo's matrix is genes x cells and alevin-fry's is cells x genes,
# so an orientation resolved the wrong way gives a plausible matrix of exactly the
# right shape carrying transposed counts; an annotation merged by position instead
# of by name gives the right columns against the wrong cells. Neither raises.
#
# tests/lib/make_annotation_fixture.py therefore writes counts that identify their
# own gene and cell, and hands the merge a CellSweep object whose barcodes are
# shuffled and whose genes are a subset. The cases below assert on those values
# rather than on shapes or exit status.
#
# No sequencing data, no mapper and no cluster are needed.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

FILTER_TOOL="$PROJECT_ROOT/bin/filter_doublets.py"
H5AD_TOOL="$PROJECT_ROOT/bin/mtx_to_h5ad.py"
TENX_TOOL="$PROJECT_ROOT/bin/mtx_to_10x.py"
FIXTURE_GEN="$TESTS_DIR/lib/make_annotation_fixture.py"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh annotated_h5ad [-- OPTIONS]
       tests/checks/annotated_h5ad.sh [OPTIONS]

Validate the count-matrix triplet handling and the annotated h5ad it ends in.

Options:
  --python PATH   Python interpreter to use (default: \$BCA_PYTHON, else
                  python3, else python).
  --keep          Keep the generated fixtures and outputs for inspection.
  -h, --help      Show this message.

Cases:
  help                the tools start and print usage
  orientation         both mappers' triplets read back as the same matrix
  dtype               the on-disk precision is preserved, and the opt-in
                      narrowing to 32 bits changes no value
  roundtrip           write_triplet -> read_triplet is value-stable
  filter_barcodes     the doublet filter keeps exactly the singlets, in order
  filter_values       each kept cell's counts survive the filter unswapped
  filter_summary      the summary counts the evaluated cells, not every droplet
  filter_alevin       the filter gives the same result on either orientation
  h5ad_plain          the bare conversion carries the counts and sample_id
  h5ad_doublets       obs['doublet_status'] flags exactly the consensus calls
  h5ad_cellsweep      the denoised layer and annotations land by name, on both
                      the matching-axes fast path and the general scatter
  h5ad_alevin         either orientation gives the same annotated object
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

log_header "annotated_h5ad"

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

# matplotlib is only in the doublet-filter environment; anndata only where the
# h5ad steps actually run. Each group of cases is gated on what it needs.
HAVE_MATPLOTLIB=0
"$PYTHON" -c 'import matplotlib' >/dev/null 2>&1 && HAVE_MATPLOTLIB=1
HAVE_ANNDATA=0
"$PYTHON" -c 'import anndata' >/dev/null 2>&1 && HAVE_ANNDATA=1

WORKDIR="$BCA_TEST_LOGDIR/annotated_h5ad"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# assert NAME PYTHON_EXPR DETAIL [PATH...]
#
# PYTHON_EXPR is evaluated with the fixture's helpers and bin/mtx_io.py in scope,
# and must yield (ok, message). The fixture module is importable so the expected
# values are named by the same rule that generated them, never restated by hand.
#
# Paths are passed as arguments rather than interpolated into the expression:
# under Git Bash a Windows interpreter receives converted paths in argv, but a
# path baked into the -c payload stays in MSYS form and cannot be opened. They
# arrive as LIB, BIN, W and then A[0], A[1], ... in the order given.
assert() {
    local name="$1" expr="$2" detail="${3:-}"; shift 3 || shift $#
    local out
    if out="$("$PYTHON" -c "
import sys
LIB, BIN, W, A = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
sys.path.insert(0, LIB)
sys.path.insert(0, BIN)
import numpy as np, pandas as pd, scipy.io as sio
import mtx_io
from make_annotation_fixture import (
    barcode, gene, raw_value, denoised_value, alpha_hat, ambient_hat, is_empty, celltype,
    barcodes, genes, doublet_barcodes, kept_barcodes, raw_matrix,
    N_GENES, N_CELLS, N_EVALUATED, DOUBLET_CELLS, CS_DROPPED_GENE,
)

def read_lines(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]

def read_summary(path):
    return dict(l.split('\t', 1) for l in read_lines(path))

def dense(matrix):
    return np.asarray(matrix.todense()) if hasattr(matrix, 'todense') else np.asarray(matrix)

def read_h5ad(path):
    import anndata as ad
    return ad.read_h5ad(path)

def read_starsolo(d, **kwargs):
    'A triplet under STARsolo file names.'
    return mtx_io.read_triplet(d + '/matrix.mtx', d + '/barcodes.tsv', d + '/features.tsv', **kwargs)

def read_alevin(d, **kwargs):
    'A triplet under alevin-fry file names, holding the same matrix transposed.'
    return mtx_io.read_triplet(d + '/quants_mat.mtx', d + '/quants_mat_rows.txt', d + '/quants_mat_cols.txt', **kwargs)

def rewrite(src, out):
    'Read a triplet and write it back out, as the doublet filter does.'
    matrix, bcs, features = read_starsolo(src)
    mtx_io.write_triplet(out, matrix, bcs, features)
    return read_starsolo(out)

ok, msg = ($expr)
print(('PASS' if ok else 'FAIL') + '\t' + str(msg))
" "$TESTS_DIR/lib" "$PROJECT_ROOT/bin" "$WORKDIR" "$@" 2>&1)"; then
        local status="${out%%$'\t'*}" msg="${out#*$'\t'}"
        if [[ "$status" == "PASS" ]]; then
            record PASS "$name" "${detail:-$msg}"
        else
            record FAIL "$name" "$msg"
        fi
    else
        record FAIL "$name" "assertion crashed: $(printf '%s' "$out" | tail -3 | tr '\n' ' ')"
    fi
}

# --------------------------------------------------------------------------
# Case: help
# --------------------------------------------------------------------------

HELP_TOOLS=("$TENX_TOOL")
[[ "$HAVE_MATPLOTLIB" -eq 1 ]] && HELP_TOOLS+=("$FILTER_TOOL")
[[ "$HAVE_ANNDATA" -eq 1 ]] && HELP_TOOLS+=("$H5AD_TOOL")

HELP_OK=1
for tool in "${HELP_TOOLS[@]}"; do
    name="$(basename "$tool")"
    if ! run_logged "$BCA_TEST_LOGDIR/annotated_h5ad_help_$(slugify "$name").log" \
            "$PYTHON" "$tool" --help; then
        record FAIL "help.$name" "--help failed"
        HELP_OK=0
    fi
done
if [[ "$HELP_OK" -eq 1 ]]; then
    record PASS "help" "${#HELP_TOOLS[@]} tools print usage"
else
    finish_check
    exit 1
fi

# --------------------------------------------------------------------------
# Fixture
# --------------------------------------------------------------------------

FIXTURE="$WORKDIR/fixture"
if ! run_logged "$BCA_TEST_LOGDIR/annotated_h5ad_fixture.log" \
        "$PYTHON" "$FIXTURE_GEN" "$FIXTURE"; then
    record FAIL "fixture" "see annotated_h5ad_fixture.log"
    finish_check
    exit 1
fi
STARSOLO="$FIXTURE/starsolo"
ALEVIN="$FIXTURE/alevin"
CALLS="$FIXTURE/combined_results.tsv"
CS_H5AD="$FIXTURE/cellsweep_full.h5ad"
CS_ALIGNED="$FIXTURE/cellsweep_aligned.h5ad"
record PASS "fixture" "both orientations, consensus calls and two CellSweep objects written"

# --------------------------------------------------------------------------
# Case: orientation
#
# The one mapper-specific decision left downstream of COLLAPSE_ALEVIN_USA. Both
# triplets hold the same counts, so they have to read back identically.
# --------------------------------------------------------------------------

assert "orientation.starsolo" \
    "(np.array_equal(dense(read_starsolo(A[0])[0]), raw_matrix()),
      'genes x cells read as given')" \
    "STARsolo's genes x cells matrix is not transposed" \
    "$STARSOLO"

assert "orientation.alevin" \
    "(np.array_equal(dense(read_alevin(A[0])[0]), raw_matrix()),
      'cells x genes transposed to genes x cells')" \
    "alevin-fry's cells x genes matrix is transposed" \
    "$ALEVIN"

assert "orientation.axes" \
    "(read_starsolo(A[0])[1] == read_alevin(A[1])[1] == barcodes()
      and list(read_starsolo(A[0])[2][0]) == list(read_alevin(A[1])[2][0]) == genes(),
      'both triplets give the same axes')" \
    "one barcode and one feature axis for both mappers" \
    "$STARSOLO" "$ALEVIN"

# --------------------------------------------------------------------------
# Cases: dtype
#
# The readers hand back the matrix at the precision it was written in. Narrowing
# it to 32 bits halves the peak memory of every step that reads a raw matrix, but
# it is opt-in per call rather than the default, so the published counts are bit
# for bit the mapper's own.
# --------------------------------------------------------------------------

assert "dtype.preserved" \
    "(read_starsolo(A[0])[0].dtype == sio.mmread(A[0] + '/matrix.mtx').dtype,
      f'counts read as {read_starsolo(A[0])[0].dtype}, as written')" \
    "the reader does not change the on-disk precision" \
    "$STARSOLO"

assert "dtype.opt_in" \
    "((lambda m: (m.dtype == np.int32 and np.array_equal(dense(m), raw_matrix()),
                  f'downcast=True gives {m.dtype}, values unchanged'))(
        read_starsolo(A[0], downcast=True)[0]))" \
    "downcast=True narrows to int32 without changing a value" \
    "$STARSOLO"

assert "dtype.float" \
    "((lambda d: (mtx_io.downcast_counts(d).dtype == np.float32
                  and abs(float(mtx_io.downcast_counts(d)[0]) - 3.4285714285714284) < 1e-6,
                  'alevin-fry EM counts keep ~7 significant digits'))(
        np.array([3.4285714285714284], dtype=np.float64)))" \
    "fractional counts would narrow to float32 within 1e-6" \
    "$STARSOLO"

# --------------------------------------------------------------------------
# Case: roundtrip
# --------------------------------------------------------------------------

assert "roundtrip" \
    "((lambda r: (np.array_equal(dense(r[0]), raw_matrix()) and r[1] == barcodes()
                  and list(r[2][0]) == genes(),
                  'values and axes survive a rewrite'))(rewrite(A[0], A[1])))" \
    "a written triplet reads back unchanged" \
    "$STARSOLO" "$WORKDIR/roundtrip"

# --------------------------------------------------------------------------
# Cases: the doublet filter
#
# It now slices the matrix rather than an AnnData object, and its output is what
# CellSweep denoises and MTX_TO_H5AD publishes.
# --------------------------------------------------------------------------

if [[ "$HAVE_MATPLOTLIB" -eq 0 ]]; then
    record SKIP "filter" "matplotlib is needed for the filtering summary plot"
else
    FILTERED="$WORKDIR/filtered_starsolo"
    if run_logged "$BCA_TEST_LOGDIR/annotated_h5ad_filter.log" \
            "$PYTHON" "$FILTER_TOOL" \
                --mtx "$STARSOLO/matrix.mtx" \
                --barcodes "$STARSOLO/barcodes.tsv" \
                --features "$STARSOLO/features.tsv" \
                --combined_results "$CALLS" \
                --method AnySinglet \
                --outdir "$FILTERED" \
                --summary_txt "$WORKDIR/filter_summary.txt" \
                --image_prefix "$WORKDIR/filter_"; then
        record PASS "filter.run" "filtered triplet written"
    else
        record FAIL "filter.run" "see annotated_h5ad_filter.log"
        finish_check
        exit 1
    fi

    assert "filter_barcodes" \
        "(read_lines(A[0] + '/barcodes.tsv') == kept_barcodes(),
          f\"got {read_lines(A[0] + '/barcodes.tsv')}\")" \
        "exactly the singlets, in source order" \
        "$FILTERED"

    # A column dropped by index instead of by barcode would still give the right
    # number of cells, with the wrong counts in them.
    assert "filter_values" \
        "((lambda m, b: (all(dense(m)[g, i] == raw_value(g, barcodes().index(bc))
                             for i, bc in enumerate(b) for g in range(N_GENES)),
                         'every kept cell keeps its own counts'))(
            read_starsolo(A[0])[0], read_starsolo(A[0])[1]))" \
        "counts follow the barcode, not the column index" \
        "$FILTERED"

    assert "filter_summary" \
        "((lambda s: (int(s['Total_Cells_In_Input_Matrix']) == N_CELLS
                      and int(s['Cells_Evaluated_For_Doublets']) == N_EVALUATED
                      and int(s['Consensus_Doublets_Removed']) == len(DOUBLET_CELLS)
                      and int(s['Cells_Remaining']) == N_CELLS - len(DOUBLET_CELLS),
                      f'summary reads {s}'))(read_summary(A[0])))" \
        "every droplet counted once, evaluated cells separately" \
        "$WORKDIR/filter_summary.txt"

    # The alevin-fry matrix arrives transposed; the filter must still slice cells.
    FILTERED_AF="$WORKDIR/filtered_alevin"
    if run_logged "$BCA_TEST_LOGDIR/annotated_h5ad_filter_alevin.log" \
            "$PYTHON" "$FILTER_TOOL" \
                --mtx "$ALEVIN/quants_mat.mtx" \
                --barcodes "$ALEVIN/quants_mat_rows.txt" \
                --features "$ALEVIN/quants_mat_cols.txt" \
                --combined_results "$CALLS" \
                --method AnySinglet \
                --outdir "$FILTERED_AF" \
                --summary_txt "$WORKDIR/filter_summary_alevin.txt" \
                --image_prefix "$WORKDIR/filter_alevin_"; then
        assert "filter_alevin" \
            "((lambda s, a: (np.array_equal(dense(s[0]), dense(a[0])) and s[1] == a[1],
                             'both orientations filter to the same matrix'))(
                read_starsolo(A[0]), read_starsolo(A[1])))" \
            "one result for both quantification methods" \
            "$FILTERED" "$FILTERED_AF"
    else
        record FAIL "filter_alevin" "see annotated_h5ad_filter_alevin.log"
    fi
fi

# --------------------------------------------------------------------------
# Cases: the annotated h5ad
# --------------------------------------------------------------------------

if [[ "$HAVE_ANNDATA" -eq 0 ]]; then
    record SKIP "h5ad" "anndata is needed to read the published object"
    finish_check
    exit 0
fi

PLAIN="$WORKDIR/plain.h5ad"
if run_logged "$BCA_TEST_LOGDIR/annotated_h5ad_plain.log" \
        "$PYTHON" "$H5AD_TOOL" \
            --mtx "$STARSOLO/matrix.mtx" \
            --barcodes "$STARSOLO/barcodes.tsv" \
            --features "$STARSOLO/features.tsv" \
            --sample_id sampleA \
            --output_h5ad "$PLAIN"; then
    assert "h5ad_plain" \
        "((lambda a: (np.array_equal(dense(a.X).T, raw_matrix())
                      and list(a.obs_names) == barcodes() and list(a.var_names) == genes()
                      and set(a.obs['sample_id']) == {'sampleA'}
                      and len(a.layers) == 0 and 'doublet_status' not in a.obs,
                      f'{a.n_obs} x {a.n_vars}, layers {list(a.layers)}'))(read_h5ad(A[0])))" \
        "cells x genes, no annotations without inputs" \
        "$PLAIN"
else
    record FAIL "h5ad_plain" "see annotated_h5ad_plain.log"
fi

ANNOTATED="$WORKDIR/annotated.h5ad"
if run_logged "$BCA_TEST_LOGDIR/annotated_h5ad_merge.log" \
        "$PYTHON" "$H5AD_TOOL" \
            --mtx "$STARSOLO/matrix.mtx" \
            --barcodes "$STARSOLO/barcodes.tsv" \
            --features "$STARSOLO/features.tsv" \
            --sample_id sampleA \
            --output_h5ad "$ANNOTATED" \
            --doublet_results "$CALLS" \
            --doublet_method AnySinglet \
            --cellsweep_h5ad "$CS_H5AD"; then
    record PASS "h5ad_merge.run" "annotated object written"

    # Doublets are called on the cell-called matrix only; every other barcode of a raw
    # matrix was never evaluated and has to come back 'singlet', not missing.
    assert "h5ad_doublets" \
        "((lambda a: ([bc for bc in a.obs_names if a.obs['doublet_status'][bc] == 'doublet'] == doublet_barcodes()
                      and a.obs['doublet_status'].isna().sum() == 0,
                      f\"flagged {[bc for bc in a.obs_names if a.obs['doublet_status'][bc] == 'doublet']}\"))(
            read_h5ad(A[0])))" \
        "exactly the consensus calls, everything else singlet" \
        "$ANNOTATED"

    # The CellSweep object's barcodes are shuffled and one gene is missing, so a
    # positional merge lands the right values on the wrong cells.
    assert "h5ad_cellsweep.layer" \
        "((lambda L: (all(L[i, g] == (0 if g == CS_DROPPED_GENE else denoised_value(g, i))
                          for i in range(N_CELLS) for g in range(N_GENES)),
                      'denoised counts placed by name'))(
            dense(read_h5ad(A[0]).layers['cellsweep'])))" \
        "every denoised count on its own cell and gene" \
        "$ANNOTATED"

    assert "h5ad_cellsweep.obs" \
        "((lambda a: (all(a.obs['alpha_hat'][barcode(c)] == alpha_hat(c) for c in range(N_CELLS))
                      and all(bool(a.obs['is_empty'][barcode(c)]) == is_empty(c) for c in range(N_CELLS))
                      and all(str(a.obs['celltype'][barcode(c)]) == celltype(c) for c in range(N_CELLS)),
                      'obs annotations reindexed by barcode'))(read_h5ad(A[0])))" \
        "alpha_hat, is_empty and celltype follow the barcode" \
        "$ANNOTATED"

    assert "h5ad_cellsweep.var" \
        "((lambda a: (all((np.isnan(a.var['ambient_hat'][gene(g)]) if g == CS_DROPPED_GENE
                           else a.var['ambient_hat'][gene(g)] == ambient_hat(g)) for g in range(N_GENES)),
                      'var annotations reindexed by gene'))(read_h5ad(A[0])))" \
        "ambient_hat follows the gene; the absent one is NaN" \
        "$ANNOTATED"

    # The raw counts are still the deliverable: merging must not overwrite X.
    assert "h5ad_cellsweep.counts" \
        "((lambda a: (np.array_equal(dense(a.X).T, raw_matrix()),
                      'X still holds the raw counts'))(read_h5ad(A[0])))" \
        "the denoised counts go in a layer, not over X" \
        "$ANNOTATED"
else
    record FAIL "h5ad_merge.run" "see annotated_h5ad_merge.log"
fi

# CellSweep reads the same triplet, so its axes normally match exactly and the merge
# takes its matrix as it stands rather than scattering it. Same answer, one copy less.
ALIGNED_OUT="$WORKDIR/annotated_aligned.h5ad"
if run_logged "$BCA_TEST_LOGDIR/annotated_h5ad_merge_aligned.log" \
        "$PYTHON" "$H5AD_TOOL" \
            --mtx "$STARSOLO/matrix.mtx" \
            --barcodes "$STARSOLO/barcodes.tsv" \
            --features "$STARSOLO/features.tsv" \
            --sample_id sampleA \
            --output_h5ad "$ALIGNED_OUT" \
            --cellsweep_h5ad "$CS_ALIGNED"; then
    assert "h5ad_cellsweep.aligned" \
        "((lambda L: (all(L[i, g] == denoised_value(g, i)
                          for i in range(N_CELLS) for g in range(N_GENES)),
                      'matching axes taken as they are'))(
            dense(read_h5ad(A[0]).layers['cellsweep'])))" \
        "the fast path places the same counts as the scatter" \
        "$ALIGNED_OUT"
else
    record FAIL "h5ad_cellsweep.aligned" "see annotated_h5ad_merge_aligned.log"
fi

ANNOTATED_AF="$WORKDIR/annotated_alevin.h5ad"
if run_logged "$BCA_TEST_LOGDIR/annotated_h5ad_merge_alevin.log" \
        "$PYTHON" "$H5AD_TOOL" \
            --mtx "$ALEVIN/quants_mat.mtx" \
            --barcodes "$ALEVIN/quants_mat_rows.txt" \
            --features "$ALEVIN/quants_mat_cols.txt" \
            --sample_id sampleA \
            --output_h5ad "$ANNOTATED_AF" \
            --doublet_results "$CALLS" \
            --doublet_method AnySinglet \
            --cellsweep_h5ad "$CS_H5AD"; then
    assert "h5ad_alevin" \
        "((lambda s, a: (np.array_equal(dense(s.X), dense(a.X))
                         and np.array_equal(dense(s.layers['cellsweep']), dense(a.layers['cellsweep']))
                         and list(s.obs_names) == list(a.obs_names)
                         and list(s.obs['doublet_status']) == list(a.obs['doublet_status']),
                         'both mappers give the same annotated object'))(
            read_h5ad(A[0]), read_h5ad(A[1])))" \
        "one published object for both quantification methods" \
        "$ANNOTATED" "$ANNOTATED_AF"
else
    record FAIL "h5ad_alevin" "see annotated_h5ad_merge_alevin.log"
fi

finish_check
