#!/usr/bin/env bash
# description: Check that the cell-calling estimators find a planted cell/ambient boundary
#
# Cell calling has no loud failure mode. Every cutoff is a plausible integer, the
# JSON always validates, and a threshold that keeps the entire ambient cloud looks
# exactly like one that does not. The pipeline has run this path untested since it
# was written, which is only survivable because the report plots the curve for a
# human to check.
#
# tests/lib/make_knee_fixture.py therefore plants the boundary: each fixture is a
# mixture of lognormal populations whose sizes are the answer, so the cases below
# assert against PLANTED_BOUNDARY rather than against a rank an earlier run
# happened to produce. The two_knee fixture carries a competing, deeper bend above
# the real cells, which is what separates the two selection rules; the no_knee
# fixture has no boundary at all, and the only correct behaviour there is to say so.
#
# No sequencing data, no STAR and no cluster are needed.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

NEW_TOOL="$PROJECT_ROOT/bin/secondderiv_cellcalling2.py"
OLD_TOOL="$PROJECT_ROOT/bin/secondderiv_cellcalling.py"
COMPARE_TOOL="$PROJECT_ROOT/bin/compare_cellcalling.py"
FIXTURE_GEN="$TESTS_DIR/lib/make_knee_fixture.py"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh cell_calling [-- OPTIONS]
       tests/checks/cell_calling.sh [OPTIONS]

Validate the second-derivative cell-calling estimators against planted boundaries.

Options:
  --python PATH   Python interpreter to use (default: \$BCA_PYTHON, else
                  python3, else python).
  --keep          Keep the generated fixtures and outputs for inspection.
  -h, --help      Show this message.

Cases:
  help              the tools start and print usage
  sharp_new         the new estimator finds an unambiguous boundary
  sharp_old         the current estimator finds it too, so the fixture is fair
  sharp_agree       the two estimators agree on the easy case
  json_contract     the JSON carries every key the dashboard reads
  deriv_thinned     the derivative arrays are paired and within --export-points
  cutoff_format     cutoff.txt is a bare integer with no trailing newline
  manual_cutoff     -m bypasses the search and still exports the curve
  two_knee_select   prominence + nearest-expected picks the real boundary
  no_knee_basin     a library with no boundary reports a loosely determined one
  tiny_fallback     a five-barcode curve degrades to the documented fallback
  compare_driver    compare_cellcalling.py runs all arms and writes its table
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

log_header "cell_calling"

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

if ! "$PYTHON" -c 'import numpy, scipy.signal' >/dev/null 2>&1; then
    record SKIP "python" "numpy and scipy are needed by the estimator"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/cell_calling"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# assert NAME PYTHON_EXPR DETAIL [PATH...]
#
# PYTHON_EXPR is evaluated with the fixture's helpers in scope and must yield
# (ok, message). The fixture module is importable so the planted boundary is
# named by the rule that generated it, never restated by hand.
#
# Paths are passed as arguments rather than interpolated into the expression:
# under Git Bash a Windows interpreter receives converted paths in argv, but a
# path baked into the -c payload stays in MSYS form and cannot be opened. They
# arrive as LIB, W and then A[0], A[1], ... in the order given.
assert() {
    local name="$1" expr="$2" detail="${3:-}"; shift 3 || shift $#
    local out
    if out="$("$PYTHON" -c "
import sys, json
LIB, W, A = sys.argv[1], sys.argv[2], sys.argv[3:]
sys.path.insert(0, LIB)
import numpy as np
from make_knee_fixture import (
    PLANTED_BOUNDARY, DECOY_BEND, EXPECTED_CELLS, TOLERANCE_DECADES, MIN_UMIS,
    ARCHETYPES, build, n_above, within_tolerance,
)

LEGACY_KEYS = ['logX', 'logY', 'customdata', 'derivX', 'derivY',
               'threshold_logX', 'threshold_logY', 'threshold_umi', 'status', 'message']

def entry(path):
    with open(path) as fh:
        payload = json.load(fh)
    return next(iter(payload.values()))

def raw(path):
    with open(path, 'rb') as fh:
        return fh.read()

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

# call NAME INPUT PREFIX [EXTRA ARGS...] -- run the new estimator on a fixture.
#
# A failing run is recorded rather than propagated: under `set -e` a non-zero exit
# here would abort the whole check, losing every case after it, when what is
# wanted is one FAIL and the remaining cases still reported.
call() {
    local name="$1" input="$2" prefix="$3"; shift 3
    local logfile="$BCA_TEST_LOGDIR/cell_calling_$(slugify "$prefix").log"
    if ! run_logged "$logfile" \
            "$PYTHON" "$NEW_TOOL" \
            -i "$WORKDIR/$input" -s "$name" \
            -o "$WORKDIR/${prefix}.json" -c "$WORKDIR/${prefix}_cutoff.txt" "$@"; then
        record FAIL "run_$prefix" "the estimator exited non-zero"
        tail_log "$logfile"
    fi
}

# --------------------------------------------------------------------------
# Case: help
# --------------------------------------------------------------------------

HELP_OK=1
for tool in "$NEW_TOOL" "$COMPARE_TOOL"; do
    if ! "$PYTHON" "$tool" --help >/dev/null 2>&1; then
        record FAIL "help" "$(basename "$tool") --help failed"
        HELP_OK=0
    fi
done
if [[ "$HELP_OK" -eq 1 ]]; then
    record PASS "help" "both tools print usage"
else
    finish_check
    exit 1
fi

# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

if ! run_logged "$BCA_TEST_LOGDIR/cell_calling_fixture.log" \
        "$PYTHON" "$FIXTURE_GEN" --outdir "$WORKDIR"; then
    record FAIL "fixture" "could not build the fixtures"
    tail_log "$BCA_TEST_LOGDIR/cell_calling_fixture.log"
    finish_check
    exit 1
fi
record PASS "fixture" "three archetypes plus a five-barcode curve"

# --------------------------------------------------------------------------
# Case: sharp -- an unambiguous boundary, which everything must find
# --------------------------------------------------------------------------

call sharp sharp_UMIperCellSorted.txt sharp -e 800

run_logged "$BCA_TEST_LOGDIR/cell_calling_sharp_old.log" \
    "$PYTHON" "$OLD_TOOL" \
    -i "$WORKDIR/sharp_UMIperCellSorted.txt" -s sharp \
    -o "$WORKDIR/sharp_old.json" -c "$WORKDIR/sharp_old_cutoff.txt" -e 800 || true

assert "sharp_new" \
    "within_tolerance(entry(A[0]).get('cutoff_rank'), 'sharp')[0], \
     'rank %s vs planted %s (%.3f decades)' % (entry(A[0]).get('cutoff_rank'), PLANTED_BOUNDARY['sharp'], within_tolerance(entry(A[0]).get('cutoff_rank'), 'sharp')[1])" \
    "" "$WORKDIR/sharp.json"

assert "sharp_old" \
    "within_tolerance(int(round(10 ** entry(A[0])['threshold_logX'])), 'sharp')[0], \
     'legacy rank %s vs planted %s' % (int(round(10 ** entry(A[0])['threshold_logX'])), PLANTED_BOUNDARY['sharp'])" \
    "" "$WORKDIR/sharp_old.json"

assert "sharp_agree" \
    "abs(entry(A[0])['threshold_logX'] - entry(A[1])['threshold_logX']) <= TOLERANCE_DECADES, \
     '%.3f decades between the new and legacy cutoffs' % abs(entry(A[0])['threshold_logX'] - entry(A[1])['threshold_logX'])" \
    "" "$WORKDIR/sharp.json" "$WORKDIR/sharp_old.json"

# --------------------------------------------------------------------------
# Case: the output contract the dashboard depends on
# --------------------------------------------------------------------------

assert "json_contract" \
    "(lambda e: (all(k in e for k in LEGACY_KEYS) and e['status'] in ('ok', 'warning') and len(e['logX']) == len(e['logY']) == len(e['customdata']), \
                 'missing: %s' % [k for k in LEGACY_KEYS if k not in e]))(entry(A[0]))" \
    "every key the report reads is present" "$WORKDIR/sharp.json"

assert "deriv_thinned" \
    "(lambda e: (len(e['derivX']) == len(e['derivY']) and 0 < len(e['derivX']) <= 2000 and len(e['raw_slope_derivY']) == len(e['derivX']), \
                 '%d derivative points' % len(e['derivX'])))(entry(A[0]))" \
    "" "$WORKDIR/sharp.json"

assert "cutoff_format" \
    "(lambda b: (b == str(int(b)).encode() and int(b) == entry(A[1])['threshold_umi'], \
                 'cutoff.txt is %r' % b))(raw(A[0]))" \
    "bare integer, no trailing newline" "$WORKDIR/sharp_cutoff.txt" "$WORKDIR/sharp.json"

# --------------------------------------------------------------------------
# Case: manual cutoff -- the search is skipped, the curve is not
# --------------------------------------------------------------------------

call sharp sharp_UMIperCellSorted.txt manual -m 777

assert "manual_cutoff" \
    "(lambda e, b: (e['threshold_umi'] == 777 and int(b) == 777 and len(e['logX']) > 0 and e['selection_rule'] == 'manual', \
                    'threshold_umi=%s, cutoff.txt=%s, rule=%s' % (e['threshold_umi'], b.decode(), e['selection_rule'])))(entry(A[0]), raw(A[1]))" \
    "" "$WORKDIR/manual.json" "$WORKDIR/manual_cutoff.txt"

# --------------------------------------------------------------------------
# Case: two competing bends -- the selection rule is the whole point
# --------------------------------------------------------------------------

call two_knee two_knee_UMIperCellSorted.txt two_knee -e 1000
call two_knee two_knee_UMIperCellSorted.txt two_knee_global -e 1000 --select global

# Both halves matter. Landing on the boundary is the claim; the global rule
# landing on the decoy instead is what shows the fixture can tell them apart. An
# earlier fixture put the decoy outside the expected-cell window, where both rules
# agreed and this case passed without exercising the choice at all.
assert "two_knee_select" \
    "(lambda e, g: (within_tolerance(e.get('cutoff_rank'), 'two_knee')[0] \
                    and not within_tolerance(g.get('cutoff_rank'), 'two_knee')[0], \
                    'prominent rank %s (planted %s), global rank %s (decoy %s), %s candidates' \
                    % (e.get('cutoff_rank'), PLANTED_BOUNDARY['two_knee'], g.get('cutoff_rank'), DECOY_BEND['two_knee'], e.get('n_candidates'))))(entry(A[0]), entry(A[1]))" \
    "" "$WORKDIR/two_knee.json" "$WORKDIR/two_knee_global.json"

# --------------------------------------------------------------------------
# Case: no boundary at all -- the answer is that there is no answer
# --------------------------------------------------------------------------

call no_knee no_knee_UMIperCellSorted.txt no_knee -e 2000

# Two things must hold. A curve with no boundary must be less well localised than
# one with an obvious boundary, or the resolution diagnostic means nothing; and no
# minimum on it may clear the prominence threshold, so the run reports that it
# found nothing defensible instead of asserting a cutoff. It still writes one --
# three downstream modules read the cutoff file -- but by the "shallow" route,
# under a warning.
assert "no_knee_basin" \
    "(lambda flat, sharp: ((flat.get('resolution_decades') or 0) > (sharp.get('resolution_decades') or 0) \
                           and flat.get('selection_rule') == 'shallow' and flat['status'] == 'warning', \
                           '%s decades indistinguishable on a featureless curve vs %s on a sharp one; status %s, rule %s' \
                           % (flat.get('resolution_decades'), sharp.get('resolution_decades'), flat['status'], flat.get('selection_rule'))))(entry(A[0]), entry(A[1]))" \
    "" "$WORKDIR/no_knee.json" "$WORKDIR/sharp.json"

# --------------------------------------------------------------------------
# Case: too little data -- the documented fallback, not a crash
# --------------------------------------------------------------------------

call tiny tiny_UMIperCellSorted.txt tiny -e 800

assert "tiny_fallback" \
    "(lambda e, b: (e['logX'] == [] and e['status'] == 'warning' and e['threshold_umi'] == 0 and int(b) == 0, \
                    'status=%s, threshold=%s' % (e['status'], e['threshold_umi'])))(entry(A[0]), raw(A[1]))" \
    "" "$WORKDIR/tiny.json" "$WORKDIR/tiny_cutoff.txt"

# --------------------------------------------------------------------------
# Case: the comparison driver runs every arm end to end
# --------------------------------------------------------------------------

if run_logged "$BCA_TEST_LOGDIR/cell_calling_compare.log" \
        "$PYTHON" "$COMPARE_TOOL" \
        "$WORKDIR/sharp_UMIperCellSorted.txt" "$WORKDIR/two_knee_UMIperCellSorted.txt" \
        --sample-id sharp --sample-id two_knee -e 800 \
        -o "$WORKDIR/comparison.tsv"; then
    assert "compare_driver" \
        "(lambda rows: (len(rows) == 3 and rows[0].split('\t')[0] == 'sample' and all(len(r.split('\t')) == len(rows[0].split('\t')) for r in rows), \
                        '%d rows, %d columns' % (len(rows) - 1, len(rows[0].split('\t')))))(open(A[0]).read().rstrip('\n').split('\n'))" \
        "" "$WORKDIR/comparison.tsv"
else
    record FAIL "compare_driver" "compare_cellcalling.py exited non-zero"
    tail_log "$BCA_TEST_LOGDIR/cell_calling_compare.log"
fi

finish_check
