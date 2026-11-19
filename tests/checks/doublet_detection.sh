#!/usr/bin/env bash
# description: Check the ported Demuxafy doublet callers, their consensus and the one-caller fallback
#
# SCRUBLET, SCDBLFINDER and COMBINE_DOUBLET_RESULTS run the scripts ported from
# Demuxafy (bin/demuxafy_scrublet.py, bin/demuxafy_scdblfinder.R,
# bin/demuxafy_combine_results.R) through bin/run_scrublet_doublet.py. This check
# runs them on a synthetic 10x matrix with known heterotypic doublets
# (tests/lib/make_doublet_fixture.py) and asserts:
#
#   - each caller writes the columns the consensus and the h5ad read, and finds
#     mostly true doublets;
#   - Scrublet's no-threshold path ends in the bootstrap re-run, and any other
#     Scrublet failure fails the wrapper instead of reading as 0% doublets;
#   - scDblFinder is seeded, so a rerun gives the same calls;
#   - the four consensus rules, with both tools and with either one alone (the
#     fallback when a caller could not run), and the errors for no input or an
#     unknown rule;
#   - filter_doublets.py and doublet_calls.annotate_callers work on a one-tool table.
#
# Needs the two module environments (modules/local/tools/scrublet/environment.yml and
# modules/local/tools/scdblfinder/environment.yml); no mapper, no Nextflow.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

BIN="$PROJECT_ROOT/bin"
FIXTURE_GEN="$TESTS_DIR/lib/make_doublet_fixture.py"

PYTHON="${BCA_SCRUBLET_PYTHON:-}"
RSCRIPT="${BCA_SCDBLFINDER_RSCRIPT:-}"
H5AD_PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh doublet_detection [-- OPTIONS]
       tests/checks/doublet_detection.sh [OPTIONS]

Run the ported doublet callers and their consensus on a synthetic matrix.

Options:
  --python PATH       Python of the scrublet environment (default:
                      \$BCA_SCRUBLET_PYTHON, else python3).
  --rscript PATH      Rscript of the scDblFinder environment (default:
                      \$BCA_SCDBLFINDER_RSCRIPT, else Rscript).
  --h5ad-python PATH  Python with anndata, for the annotate case (default:
                      \$BCA_PYTHON, else --python; skipped without anndata).
  --keep              Keep the fixture and outputs for inspection.
  -h, --help          Show this message.

Cases:
  help                every script starts and prints usage
  scrublet            Scrublet's columns, and mostly true doublets called
  scrublet_rate       the expected rate follows --doublet_rate_per_1000, capped at 0.5
  scrublet_fallback   no threshold: exit 3, and the wrapper's bootstrap re-run
  scrublet_error      another failure fails the wrapper
  scdblfinder         scDblFinder's columns, mostly true doublets, same calls on a rerun
  combine             the four rules with both tools: intersection and union
  combine_single      with one tool, every rule gives that tool's calls
  combine_errors      no input, or an unknown rule, fails
  filter_single       filter_doublets.py on a one-tool table
  annotate_single     .uns['doublet_callers'] lists the one tool
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --python)      PYTHON="${2:?--python needs a value}"; shift 2 ;;
        --rscript)     RSCRIPT="${2:?--rscript needs a value}"; shift 2 ;;
        --h5ad-python) H5AD_PYTHON="${2:?--h5ad-python needs a value}"; shift 2 ;;
        --keep)        KEEP=1; shift ;;
        -h|--help)     usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

# --------------------------------------------------------------------------
# Pre-flight
# --------------------------------------------------------------------------

log_header "doublet_detection"

[[ -n "$PYTHON" ]] || PYTHON="python3"
[[ -n "$RSCRIPT" ]] || RSCRIPT="Rscript"
[[ -n "$H5AD_PYTHON" ]] || H5AD_PYTHON="$PYTHON"

if ! "$PYTHON" -c 'import scrublet, skimage, matplotlib, pandas, scipy' >/dev/null 2>&1; then
    record SKIP "python" "no python with scrublet (pass the scrublet environment's python with --python)"
    finish_check
    exit 0
fi
if ! "$RSCRIPT" -e 'suppressMessages({library(scDblFinder); library(argparse); library(data.table)})' >/dev/null 2>&1; then
    record SKIP "rscript" "no Rscript with scDblFinder (pass the scDblFinder environment's Rscript with --rscript)"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/doublet_detection"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# The scripts call each other by name (run_scrublet_doublet.py -> demuxafy_scrublet.py)
# and through their shebangs, so each environment goes first on PATH, as in the
# module's container. numba and matplotlib cache into the work directory.
PY_PATH="$(dirname "$PYTHON"):$BIN:$PATH"
R_PATH="$(dirname "$RSCRIPT"):$BIN:$PATH"
export NUMBA_CACHE_DIR="$WORKDIR/.numba_cache" MPLCONFIGDIR="$WORKDIR/.mplconfig"
py() { PATH="$PY_PATH" "$@"; }
rs() { PATH="$R_PATH" "$@"; }

# assert NAME PYTHON_EXPR DETAIL [PATH...]
#
# PYTHON_EXPR yields (ok, message), with pandas as pd, the fixture's truth as TRUTH and
# the paths as A[0], A[1], ... in scope; calls(path, col) gives the barcodes a table
# calls doublets in that column.
assert() {
    local name="$1" expr="$2" detail="${3:-}"; shift 3 || shift $#
    local out
    if out="$("$PYTHON" -c "
import sys
import pandas as pd
W, A = sys.argv[1], sys.argv[2:]
TRUTH = set(open(W + '/fx/mixed_doublets.txt').read().split())
def table(path):
    return pd.read_csv(path, sep='\t', dtype={'Barcode': str})
def calls(path, col):
    t = table(path)
    return set(t.loc[t[col].astype(str).str.lower() == 'doublet', 'Barcode'])
ok, msg = ($expr)
print(('PASS' if ok else 'FAIL') + '\t' + str(msg))
" "$WORKDIR" "$@" 2>&1)"; then
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

# --------------------------------------------------------------------------
# Case: help
# --------------------------------------------------------------------------

HELP_OK=1
for tool in demuxafy_scrublet.py run_scrublet_doublet.py; do
    run_logged "$BCA_TEST_LOGDIR/doublet_detection_help_$tool.log" py "$PYTHON" "$BIN/$tool" --help \
        || { record FAIL "help.$tool" "--help failed"; HELP_OK=0; }
done
for tool in demuxafy_scdblfinder.R demuxafy_combine_results.R; do
    run_logged "$BCA_TEST_LOGDIR/doublet_detection_help_$tool.log" rs "$RSCRIPT" "$BIN/$tool" --help \
        || { record FAIL "help.$tool" "--help failed"; HELP_OK=0; }
done
if [[ "$HELP_OK" -eq 1 ]]; then
    record PASS "help" "4 scripts print usage"
else
    finish_check
    exit 1
fi

# --------------------------------------------------------------------------
# Fixture
# --------------------------------------------------------------------------

if ! run_logged "$BCA_TEST_LOGDIR/doublet_detection_fixture.log" "$PYTHON" "$FIXTURE_GEN" "$WORKDIR/fx"; then
    record FAIL "fixture" "see doublet_detection_fixture.log"
    finish_check
    exit 1
fi
MIXED="$WORKDIR/fx/mixed"
record PASS "fixture" "660 cells (60 heterotypic doublets) and a homogeneous matrix"

# --------------------------------------------------------------------------
# Case: scrublet
# --------------------------------------------------------------------------

SCR="$WORKDIR/scrublet/scrublet_results.tsv"
if run_logged "$BCA_TEST_LOGDIR/doublet_detection_scrublet.log" \
        py "$PYTHON" "$BIN/run_scrublet_doublet.py" --counts_matrix "$MIXED" \
        --outdir "$WORKDIR/scrublet" --sample_id mixed; then
    assert "scrublet" \
        "(lambda t, c: (list(t.columns) == ['Barcode', 'scrublet_DropletType', 'scrublet_Scores']
                        and len(t) == 660 and len(c) >= 30 and len(c & TRUTH) >= 0.8 * len(c),
                        f'{len(c)} called, {len(c & TRUTH)} true, columns {list(t.columns)}'))
         (table(A[0]), calls(A[0], 'scrublet_DropletType'))" \
        "" "$SCR"
else
    record FAIL "scrublet" "run_scrublet_doublet.py failed, see doublet_detection_scrublet.log"
fi

# --------------------------------------------------------------------------
# Case: scrublet_rate
# --------------------------------------------------------------------------

rate_log="$BCA_TEST_LOGDIR/doublet_detection_scrublet_rate.log"
cap_log="$BCA_TEST_LOGDIR/doublet_detection_scrublet_cap.log"
if run_logged "$rate_log" py "$PYTHON" "$BIN/demuxafy_scrublet.py" -m "$MIXED" -o "$WORKDIR/rate" \
        --doublet_rate_per_1000 0.02 \
   && run_logged "$cap_log" py "$PYTHON" "$BIN/demuxafy_scrublet.py" -m "$MIXED" -o "$WORKDIR/cap" \
        --doublet_rate_per_1000 1; then
    # 660 cells: 0.02 gives 0.0132; 1 gives 0.66, capped at 0.5
    if grep -q "Expected doublet rate: 0.0132" "$rate_log" \
       && grep -q "capped at 0.5" "$cap_log" && grep -q "Expected doublet rate: 0.5000" "$cap_log"; then
        record PASS "scrublet_rate" "660 cells x 0.02/1000 = 0.0132; x 1/1000 capped at 0.5"
    else
        record FAIL "scrublet_rate" "expected rate not in the logs (doublet_detection_scrublet_rate/cap.log)"
    fi
else
    record FAIL "scrublet_rate" "demuxafy_scrublet.py failed"
fi

# --------------------------------------------------------------------------
# Case: scrublet_fallback
#
# Scrublet sets no threshold when its simulated scores have no second mode. Which
# real matrices do that is hard to pin down, so a shim named demuxafy_scrublet.py,
# first on PATH, makes skimage's threshold_minimum raise as it does then.
# --------------------------------------------------------------------------

SHIM="$WORKDIR/shim"
mkdir -p "$SHIM"
cat > "$SHIM/demuxafy_scrublet.py" <<EOF
#!$PYTHON
import runpy, sys
import skimage.filters
def no_split(*args, **kwargs):
    raise RuntimeError("forced: no bimodal split")
skimage.filters.threshold_minimum = no_split
sys.argv[0] = "$BIN/demuxafy_scrublet.py"
runpy.run_path(sys.argv[0], run_name="__main__")
EOF
chmod +x "$SHIM/demuxafy_scrublet.py"

set +e
PATH="$SHIM:$PY_PATH" "$SHIM/demuxafy_scrublet.py" -m "$MIXED" -o "$WORKDIR/nothr_direct" \
    > "$BCA_TEST_LOGDIR/doublet_detection_scrublet_nothr.log" 2>&1
direct_exit=$?
set -e
fallback_log="$BCA_TEST_LOGDIR/doublet_detection_scrublet_fallback.log"
if [[ "$direct_exit" -ne 3 ]]; then
    record FAIL "scrublet_fallback" "demuxafy_scrublet.py exited $direct_exit without a threshold, expected 3"
elif run_logged "$fallback_log" env PATH="$SHIM:$PY_PATH" "$PYTHON" "$BIN/run_scrublet_doublet.py" \
        --counts_matrix "$MIXED" --outdir "$WORKDIR/fallback" --sample_id mixed \
   && grep -q "no bimodal" "$fallback_log"; then
    assert "scrublet_fallback" \
        "(lambda t: (len(t) == 660 and (t['scrublet_DropletType'] == 'singlet').all(),
                     f'exit 3, then the bootstrap re-run: {len(t)} cells, 0% doublets reported as is'))
         (table(A[0]))" \
        "" "$WORKDIR/fallback/scrublet_results.tsv"
else
    record FAIL "scrublet_fallback" "the wrapper did not recover, see doublet_detection_scrublet_fallback.log"
fi

# --------------------------------------------------------------------------
# Case: scrublet_error
# --------------------------------------------------------------------------

if py "$PYTHON" "$BIN/run_scrublet_doublet.py" --counts_matrix "$WORKDIR/no_such_matrix" \
        --outdir "$WORKDIR/error" --sample_id missing > "$BCA_TEST_LOGDIR/doublet_detection_scrublet_error.log" 2>&1; then
    record FAIL "scrublet_error" "the wrapper succeeded on a missing matrix"
elif [[ -f "$WORKDIR/error/scrublet_results.tsv" ]]; then
    record FAIL "scrublet_error" "a missing matrix still produced results"
else
    record PASS "scrublet_error" "a missing matrix fails the wrapper, without a bootstrap re-run"
fi

# --------------------------------------------------------------------------
# Case: scdblfinder
# --------------------------------------------------------------------------

SCD="$WORKDIR/scdblfinder/scDblFinder_doublets_singlets.tsv"
if run_logged "$BCA_TEST_LOGDIR/doublet_detection_scdblfinder.log" \
        rs "$RSCRIPT" "$BIN/demuxafy_scdblfinder.R" -o "$WORKDIR/scdblfinder" -t "$MIXED" \
   && run_logged "$BCA_TEST_LOGDIR/doublet_detection_scdblfinder_rerun.log" \
        rs "$RSCRIPT" "$BIN/demuxafy_scdblfinder.R" -o "$WORKDIR/scdblfinder_rerun" -t "$MIXED"; then
    if grep -q "training iteration failed" "$BCA_TEST_LOGDIR/doublet_detection_scdblfinder.log"; then
        record FAIL "scdblfinder" "scDblFinder's training failed (xgboost API mismatch?), see doublet_detection_scdblfinder.log"
    else
        assert "scdblfinder" \
            "(lambda t, c, rerun: (list(t.columns) == ['Barcode', 'scDblFinder_DropletType', 'scDblFinder_Score']
                                   and len(t) == 660 and len(c) >= 10 and len(c & TRUTH) >= 0.8 * len(c)
                                   and t.equals(rerun),
                                   f'{len(c)} called, {len(c & TRUTH)} true; rerun identical: {t.equals(rerun)}'))
             (table(A[0]), calls(A[0], 'scDblFinder_DropletType'), table(A[1]))" \
            "" "$SCD" "$WORKDIR/scdblfinder_rerun/scDblFinder_doublets_singlets.tsv"
    fi
else
    record FAIL "scdblfinder" "demuxafy_scdblfinder.R failed, see doublet_detection_scdblfinder.log"
fi

# --------------------------------------------------------------------------
# Case: combine
# --------------------------------------------------------------------------

METHODS=(AnySinglet AnyDoublet MajoritySinglet AtLeastHalfSinglet)

# combine LABEL METHOD ARGS... -- writes $WORKDIR/combine_LABEL_METHOD/c_w_combined_assignments.tsv
combine() {
    local label="$1" method="$2"; shift 2
    rs "$RSCRIPT" "$BIN/demuxafy_combine_results.R" -o "$WORKDIR/combine_${label}_${method}/c.tsv" "$@" -m "$method" \
        > "$BCA_TEST_LOGDIR/doublet_detection_combine_${label}_${method}.log" 2>&1
}

if [[ -f "$SCR" && -f "$SCD" ]]; then
    ok=1
    for m in "${METHODS[@]}"; do combine both "$m" -r "$SCR" -n "$SCD" || ok=0; done
    if [[ "$ok" -eq 1 ]]; then
        assert "combine" \
            "(lambda s, d, t, rules: (
                rules['AnySinglet'] == s & d and rules['AnyDoublet'] == s | d
                and rules['MajoritySinglet'] == rules['AnyDoublet'] and rules['AtLeastHalfSinglet'] == rules['AnySinglet']
                and list(t.columns) == ['Barcode', 'scDblFinder_DropletType', 'scDblFinder_Score',
                                        'scrublet_DropletType', 'scrublet_Scores', 'AnySinglet_DropletType'],
                f'AnySinglet = intersection ({len(s & d)}), AnyDoublet = union ({len(s | d)}); '
                f'Majority = AnyDoublet, AtLeastHalf = AnySinglet; Demuxafy column order'))
             (calls(A[0], 'scrublet_DropletType'), calls(A[1], 'scDblFinder_DropletType'), table(A[2]),
              {m: calls(A[2].replace('AnySinglet', m), m + '_DropletType')
               for m in ['AnySinglet', 'AnyDoublet', 'MajoritySinglet', 'AtLeastHalfSinglet']})" \
            "" "$SCR" "$SCD" "$WORKDIR/combine_both_AnySinglet/c_w_combined_assignments.tsv"
    else
        record FAIL "combine" "demuxafy_combine_results.R failed, see doublet_detection_combine_both_*.log"
    fi

    # ----------------------------------------------------------------------
    # Case: combine_single
    # ----------------------------------------------------------------------

    ok=1
    for m in "${METHODS[@]}"; do
        combine scr "$m" -r "$SCR" || ok=0
        combine scd "$m" -n "$SCD" || ok=0
    done
    if [[ "$ok" -eq 1 ]]; then
        assert "combine_single" \
            "(lambda s, d, ok: (ok, 'with -r or -n alone, all four rules give that tool\'s own calls, and only its columns'))
             (calls(A[0], 'scrublet_DropletType'), calls(A[1], 'scDblFinder_DropletType'),
              all(calls(A[2] + '/combine_scr_' + m + '/c_w_combined_assignments.tsv', m + '_DropletType') == calls(A[0], 'scrublet_DropletType')
                  and calls(A[2] + '/combine_scd_' + m + '/c_w_combined_assignments.tsv', m + '_DropletType') == calls(A[1], 'scDblFinder_DropletType')
                  and 'scDblFinder_DropletType' not in table(A[2] + '/combine_scr_' + m + '/c_w_combined_assignments.tsv').columns
                  and 'scrublet_DropletType' not in table(A[2] + '/combine_scd_' + m + '/c_w_combined_assignments.tsv').columns
                  for m in ['AnySinglet', 'AnyDoublet', 'MajoritySinglet', 'AtLeastHalfSinglet']))" \
            "" "$SCR" "$SCD" "$WORKDIR"
    else
        record FAIL "combine_single" "demuxafy_combine_results.R failed on one input, see doublet_detection_combine_sc*.log"
    fi
else
    record SKIP "combine" "needs both callers' results"
    record SKIP "combine_single" "needs both callers' results"
fi

# --------------------------------------------------------------------------
# Case: combine_errors
# --------------------------------------------------------------------------

if combine none AnySinglet; then
    record FAIL "combine_errors" "no input succeeded"
elif [[ -f "$SCR" ]] && combine bad NoSuchRule -r "$SCR"; then
    record FAIL "combine_errors" "an unknown rule succeeded"
else
    record PASS "combine_errors" "no input, and an unknown rule, exit non-zero"
fi

# --------------------------------------------------------------------------
# Case: filter_single
#
# The doublet filter's plot used to read both tools' columns; with one tool's table
# it must still filter and plot.
# --------------------------------------------------------------------------

ONE="$WORKDIR/combine_scd_AnySinglet/c_w_combined_assignments.tsv"
TRIPLET="$WORKDIR/triplet"
if [[ ! -f "$ONE" ]]; then
    record SKIP "filter_single" "no one-tool table"
elif ! "$PYTHON" -c "import matplotlib" >/dev/null 2>&1; then
    record SKIP "filter_single" "matplotlib missing"
else
    "$PYTHON" -W ignore - "$MIXED" "$TRIPLET" <<'EOF'
import gzip, os, sys
import scipy.io
src, out = sys.argv[1:]
os.makedirs(out, exist_ok=True)
with gzip.open(src + "/matrix.mtx.gz") as f:
    scipy.io.mmwrite(out + "/matrix.mtx", scipy.io.mmread(f))
for name in ("barcodes", "features"):
    with gzip.open(f"{src}/{name}.tsv.gz", "rt") as f, open(f"{out}/{name}.tsv", "w") as g:
        g.write(f.read())
EOF
    if run_logged "$BCA_TEST_LOGDIR/doublet_detection_filter_single.log" \
            env PYTHONPATH="$BIN" MPLBACKEND=Agg "$PYTHON" "$BIN/filter_doublets.py" \
            --mtx "$TRIPLET/matrix.mtx" --barcodes "$TRIPLET/barcodes.tsv" --features "$TRIPLET/features.tsv" \
            --combined_results "$ONE" --method AnySinglet --outdir "$WORKDIR/filtered" \
            --summary_txt "$WORKDIR/filter_summary.txt" --image_prefix "$WORKDIR/single_"; then
        assert "filter_single" \
            "(lambda kept, d: (len(kept) == 660 - len(d) and not (set(kept) & d),
                               f'{len(d)} scDblFinder doublets removed, {len(kept)} kept, plot written'))
             ([l.strip() for l in open(A[0]) if l.strip()], calls(A[1], 'AnySinglet_DropletType'))" \
            "" "$WORKDIR/filtered/barcodes.tsv" "$ONE"
        [[ -f "$WORKDIR/single_doublet_filtering_summary.png" ]] \
            || record FAIL "filter_single.plot" "no summary plot"
    else
        record FAIL "filter_single" "filter_doublets.py failed on a one-tool table, see doublet_detection_filter_single.log"
    fi
fi

# --------------------------------------------------------------------------
# Case: annotate_single
# --------------------------------------------------------------------------

if [[ ! -f "$ONE" ]]; then
    record SKIP "annotate_single" "no one-tool table"
elif ! "$H5AD_PYTHON" -c "import anndata" >/dev/null 2>&1; then
    record SKIP "annotate_single" "no python with anndata (pass one with --h5ad-python)"
else
    if out="$("$H5AD_PYTHON" - "$BIN" "$ONE" 2>&1 <<'EOF'
import sys
sys.path.insert(0, sys.argv[1])
import anndata as ad, numpy as np, pandas as pd
import doublet_calls
combined = doublet_calls.read_combined_results(sys.argv[2])
a = ad.AnnData(np.zeros((len(combined), 1)), obs=pd.DataFrame(index=combined["Barcode"].to_numpy()))
a = doublet_calls.annotate_callers(a, combined, "AnySinglet")
ok = list(a.uns["doublet_callers"]) == ["scdblfinder"] and "scdblfinder_call" in a.obs and "scrublet_call" not in a.obs
print(("PASS" if ok else "FAIL") + f"\tuns['doublet_callers'] = {list(a.uns['doublet_callers'])}, obs = {list(a.obs.columns)}")
EOF
)"; then
        last="${out##*$'\n'}"
        if [[ "$last" == PASS* ]]; then
            record PASS "annotate_single" "${last#*$'\t'}"
        else
            record FAIL "annotate_single" "${last#*$'\t'}"
        fi
    else
        record FAIL "annotate_single" "crashed: $(printf '%s' "$out" | tail -3 | tr '\n' ' ')"
    fi
fi

finish_check
