#!/usr/bin/env bash
# description: Check that the intronic / RNA-velocity matrices are subsetted and packed correctly
#
# The velocity outputs are three matrices that share one barcode and one feature
# axis, and every failure mode here is silent. Subsetting them on a UMI cutoff of
# their own instead of the GeneFull_Ex50pAS cell call would give a plausible
# matrix describing the wrong cells; transposing one and not the others, or
# mapping alevin-fry's -S block onto the unspliced layer, gives an object of
# exactly the right shape carrying the wrong numbers. None of that raises.
#
# tests/lib/make_velocyto_fixture.py therefore writes counts that identify their
# own layer, gene and cell, and the cases below assert on those values rather
# than on shapes or exit status.
#
# No sequencing data, no STAR and no cluster are needed.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

SUBSET_TOOL="$PROJECT_ROOT/bin/subset_matrices_to_cells.py"
H5AD_TOOL="$PROJECT_ROOT/bin/velocity_matrices_to_h5ad.py"
COLLAPSE_TOOL="$PROJECT_ROOT/bin/collapse_alevin_usa.py"
FIXTURE_GEN="$TESTS_DIR/lib/make_velocyto_fixture.py"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh velocity_matrix [-- OPTIONS]
       tests/checks/velocity_matrix.sh [OPTIONS]

Validate the intronic / RNA-velocity matrix outputs.

Options:
  --python PATH   Python interpreter to use (default: \$BCA_PYTHON, else
                  python3, else python).
  --keep          Keep the generated fixtures and outputs for inspection.
  -h, --help      Show this message.

Cases:
  help              the tools start and print usage
  subset_barcodes   the subset keeps exactly the called cells, in order
  subset_values     each layer's counts survive the subset unswapped
  h5ad_layers       the AnnData object carries the three layers and their sum
  alevin_unspliced  --counts U returns exactly the -U block
  absent            a Solo.out without Velocyto/ fails loudly, not silently
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

log_header "velocity_matrix"

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

WORKDIR="$BCA_TEST_LOGDIR/velocity_matrix"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# assert NAME PYTHON_EXPR DETAIL [PATH...]
#
# PYTHON_EXPR is evaluated with the fixture's helpers in scope and must yield
# (ok, message). The fixture module is importable so the expected counts are
# named by the same rule that generated them, never restated by hand.
#
# Paths are passed as arguments rather than interpolated into the expression:
# under Git Bash a Windows interpreter receives converted paths in argv, but a
# path baked into the -c payload stays in MSYS form and cannot be opened. They
# arrive as LIB, W and then A[0], A[1], ... in the order given.
assert() {
    local name="$1" expr="$2" detail="${3:-}"; shift 3 || shift $#
    local out
    if out="$("$PYTHON" -c "
import sys
LIB, W, A = sys.argv[1], sys.argv[2], sys.argv[3:]
sys.path.insert(0, LIB)
import numpy as np, pandas as pd, scipy.io as sio
from make_velocyto_fixture import layer_value, barcode, gene, N_GENES, N_CELLS, N_CALLED_CELLS

def read_lines(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]

def read_mtx(path):
    return sio.mmread(path).toarray()

ok, msg = ($expr)
print(('PASS' if ok else 'FAIL') + '\t' + str(msg))
" "$TESTS_DIR/lib" "$WORKDIR" "$@" 2>&1)"; then
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

# anndata is only present where the velocity h5ad step actually runs (the scanpy
# environment), so every case that needs it is gated on it being importable here.
HAVE_ANNDATA=0
"$PYTHON" -c 'import anndata' >/dev/null 2>&1 && HAVE_ANNDATA=1

HELP_TOOLS=("$SUBSET_TOOL" "$COLLAPSE_TOOL")
[[ "$HAVE_ANNDATA" -eq 1 ]] && HELP_TOOLS+=("$H5AD_TOOL")

HELP_OK=1
for tool in "${HELP_TOOLS[@]}"; do
    name="$(basename "$tool")"
    if ! run_logged "$BCA_TEST_LOGDIR/velocity_matrix_help_$(slugify "$name").log" \
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
if ! run_logged "$BCA_TEST_LOGDIR/velocity_matrix_fixture.log" \
        "$PYTHON" "$FIXTURE_GEN" "$FIXTURE"; then
    record FAIL "fixture" "see velocity_matrix_fixture.log"
    finish_check
    exit 1
fi
VELO_RAW="$FIXTURE/sampleA_starsolo_Solo.out/Velocyto/raw"
CALLED="$FIXTURE/called_barcodes.tsv"
ALEVIN="$FIXTURE/alevin"
record PASS "fixture" "synthetic Velocyto, GeneFull and USA matrices written"

# --------------------------------------------------------------------------
# Case: subset_barcodes
#
# The cell set comes from the GeneFull_Ex50pAS call, not from a cutoff recomputed
# on the velocity matrices. Getting this wrong is the whole reason the subset
# script takes a barcode list rather than a threshold.
# --------------------------------------------------------------------------

SUBSET_OUT="$WORKDIR/velocyto_filtered"
if run_logged "$BCA_TEST_LOGDIR/velocity_matrix_subset.log" \
        "$PYTHON" "$SUBSET_TOOL" --dir "$VELO_RAW" --barcodes "$CALLED" --outdir "$SUBSET_OUT"; then
    record PASS "subset_barcodes.run" "subset written"
else
    record FAIL "subset_barcodes.run" "see velocity_matrix_subset.log"
    finish_check
    exit 1
fi

assert "subset_barcodes.identity" \
    "(read_lines(A[0] + '/barcodes.tsv') == read_lines(A[1]),
      f\"got {read_lines(A[0] + '/barcodes.tsv')}\")" \
    "exactly the called barcodes, in the called order" \
    "$SUBSET_OUT" "$CALLED"

assert "subset_barcodes.features" \
    "(read_lines(A[0] + '/features.tsv') == [gene(i) for i in range(N_GENES)],
      'feature axis preserved')" \
    "the feature axis is untouched" \
    "$SUBSET_OUT"

# All three matrices must be cut to the same columns, or the layers stop lining up.
assert "subset_barcodes.shapes" \
    "(all(read_mtx(A[0] + f'/{l}.mtx').shape == (N_GENES, N_CALLED_CELLS)
          for l in ('spliced', 'unspliced', 'ambiguous')),
      {l: read_mtx(A[0] + f'/{l}.mtx').shape for l in ('spliced', 'unspliced', 'ambiguous')})" \
    "every layer cut to the same cells" \
    "$SUBSET_OUT"

# --------------------------------------------------------------------------
# Case: subset_values
#
# Shapes alone cannot tell a correct subset from one that kept the wrong columns
# or swapped two layers, so the counts themselves are checked against the rule
# that generated them.
# --------------------------------------------------------------------------

assert "subset_values.counts" \
    "(all(read_mtx(A[0] + f'/{l}.mtx')[g][c] == layer_value(l, g, c)
          for l in ('spliced', 'unspliced', 'ambiguous')
          for g in range(N_GENES) for c in range(N_CALLED_CELLS)),
      'a retained count does not match the value it was generated with')" \
    "each layer's counts survive unswapped" \
    "$SUBSET_OUT"

# --------------------------------------------------------------------------
# Case: h5ad_layers
# --------------------------------------------------------------------------

if [[ "$HAVE_ANNDATA" -eq 0 ]]; then
    record SKIP "h5ad_layers" "anndata is not installed"
else
    STAR_H5AD="$WORKDIR/starsolo_velocity.h5ad"
    if run_logged "$BCA_TEST_LOGDIR/velocity_matrix_h5ad_starsolo.log" \
            "$PYTHON" "$H5AD_TOOL" --starsolo-dir "$SUBSET_OUT" --out "$STAR_H5AD" --sample-id sampleA; then
        record PASS "h5ad_layers.starsolo_run" "AnnData written from the STARsolo layout"
    else
        record FAIL "h5ad_layers.starsolo_run" "see velocity_matrix_h5ad_starsolo.log"
    fi

    ALEVIN_H5AD="$WORKDIR/alevin_velocity.h5ad"
    if run_logged "$BCA_TEST_LOGDIR/velocity_matrix_h5ad_alevin.log" \
            "$PYTHON" "$H5AD_TOOL" --alevin-dir "$ALEVIN" --out "$ALEVIN_H5AD" --sample-id sampleA; then
        record PASS "h5ad_layers.alevin_run" "AnnData written from the alevin-fry USA layout"
    else
        record FAIL "h5ad_layers.alevin_run" "see velocity_matrix_h5ad_alevin.log"
    fi

    # STARsolo writes genes x cells and alevin-fry cells x genes; both must come
    # out cells x genes, or every downstream velocity tool reads the object sideways.
    assert "h5ad_layers.orientation" \
        "(__import__('anndata').read_h5ad(A[0]).shape == (N_CALLED_CELLS, N_GENES)
          and __import__('anndata').read_h5ad(A[1]).shape == (N_CELLS, N_GENES),
          f\"starsolo={__import__('anndata').read_h5ad(A[0]).shape} \"
          f\"alevin={__import__('anndata').read_h5ad(A[1]).shape}\")" \
        "both mappers' objects come out cells x genes" \
        "$STAR_H5AD" "$ALEVIN_H5AD"

    # The layer a count lands in is the entire point; a transposed read or a
    # mis-split USA block would still produce three layers of the right shape.
    assert "h5ad_layers.values" \
        "(all(
            __import__('anndata').read_h5ad(path).layers[l].toarray()[c][g] == layer_value(l, g, c)
            for path, ncells in ((A[0], N_CALLED_CELLS), (A[1], N_CELLS))
            for l in ('spliced', 'unspliced', 'ambiguous')
            for g in range(N_GENES) for c in range(ncells)),
          'a layer holds a count belonging to another layer, gene or cell')" \
        "each count lands in the right layer for both mappers" \
        "$STAR_H5AD" "$ALEVIN_H5AD"

    assert "h5ad_layers.x_is_sum" \
        "(np.array_equal(
            __import__('anndata').read_h5ad(A[0]).X.toarray(),
            sum(__import__('anndata').read_h5ad(A[0]).layers[l].toarray()
                for l in ('spliced', 'unspliced', 'ambiguous'))),
          'X is not the sum of the three layers')" \
        "X is the sum of the three layers" \
        "$STAR_H5AD"
fi

# --------------------------------------------------------------------------
# Case: alevin_unspliced
#
# alevin-fry's counterpart of the Velocyto unspliced matrix. The -S block is
# written first, so a collapse that ignored the suffix would return it instead
# and every value would still look like a plausible count.
# --------------------------------------------------------------------------

UNSPLICED_OUT="$WORKDIR/alevin_unspliced"
if run_logged "$BCA_TEST_LOGDIR/velocity_matrix_alevin_u.log" \
        "$PYTHON" "$COLLAPSE_TOOL" --dir "$ALEVIN" --outdir "$UNSPLICED_OUT" --counts U; then
    record PASS "alevin_unspliced.run" "--counts U accepted"
else
    record FAIL "alevin_unspliced.run" "see velocity_matrix_alevin_u.log"
fi

assert "alevin_unspliced.values" \
    "(all(read_mtx(A[0] + '/quants_mat.mtx')[c][g] == layer_value('unspliced', g, c)
          for g in range(N_GENES) for c in range(N_CELLS)),
      'the collapsed matrix does not hold the -U block')" \
    "exactly the unspliced block, not the spliced one" \
    "$UNSPLICED_OUT"

assert "alevin_unspliced.genes" \
    "(read_lines(A[0] + '/quants_mat_cols.txt') == [gene(i) for i in range(N_GENES)],
      'gene names are not the plain, unsuffixed IDs')" \
    "USA suffixes stripped from the gene names" \
    "$UNSPLICED_OUT"

# --------------------------------------------------------------------------
# Case: absent
#
# Most runs leave perform_velocity off, so no Velocyto directory exists. The
# subset must then refuse outright rather than write a partial or empty result
# that would be published as a velocity matrix.
# --------------------------------------------------------------------------

MISSING="$FIXTURE/sampleB_starsolo_Solo.out/Velocyto/raw"
if run_logged "$BCA_TEST_LOGDIR/velocity_matrix_absent.log" \
        "$PYTHON" "$SUBSET_TOOL" --dir "$MISSING" --barcodes "$CALLED" --outdir "$WORKDIR/absent_out"; then
    record FAIL "absent" "the subset succeeded against a missing Velocyto directory"
elif [[ -d "$WORKDIR/absent_out" ]]; then
    record FAIL "absent" "the subset failed but still created an output directory"
else
    record PASS "absent" "a missing Velocyto directory fails without writing output"
fi

finish_check
