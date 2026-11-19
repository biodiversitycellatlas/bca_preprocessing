#!/usr/bin/env bash
# description: Check that the cell caller finds a planted cell/ambient boundary
#
# Cell calling has no loud failure mode. Every cutoff is a plausible integer, the
# JSON always validates, and a threshold that keeps the entire ambient cloud looks
# exactly like one that does not. Without these checks the only safeguard is the
# report's barcode-rank plot and a human looking at it.
#
# The fixtures therefore plant the boundary, and the cases assert against the rule
# that generated it rather than against a rank an earlier run happened to produce:
#
#   tests/lib/make_knee_fixture.py     curves only (UMIperCellSorted.txt), a mixture
#                                      of lognormal populations whose sizes are the
#                                      answer. Run without --raw-dir, they test the
#                                      curve rule: a clean boundary, one below a
#                                      deeper decoy bend, an ambient shoulder, no
#                                      boundary at all, a curve too short to read.
#   tests/lib/make_ambient_fixture.py  STARsolo raw/ directories with planted
#                                      expression, for the ambient test that --raw-dir
#                                      enables: a subpopulation bend above the cells,
#                                      a debris shoulder below them, cells that thin
#                                      out with no bend, and split-pool libraries with
#                                      no empty plateau.
#
# The matrix fixtures are only worth anything if the curve alone gets them wrong, so
# the *_fair cases run the same script on them without --raw-dir and require it to
# fail. No sequencing data, no STAR and no cluster are needed.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

TOOL="$PROJECT_ROOT/bin/secondderiv_cellcalling.py"
ALEVIN_TOOL="$PROJECT_ROOT/bin/secondderiv_alevin.py"
FIXTURE_GEN="$TESTS_DIR/lib/make_knee_fixture.py"
AMBIENT_FIXTURE_GEN="$TESTS_DIR/lib/make_ambient_fixture.py"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh cell_calling [-- OPTIONS]
       tests/checks/cell_calling.sh [OPTIONS]

Validate secondderiv_cellcalling.py against planted cell/ambient boundaries.

Options:
  --python PATH   Python interpreter to use (default: \$BCA_PYTHON, else
                  python3, else python).
  --keep          Keep the generated fixtures and outputs for inspection.
  -h, --help      Show this message.

Cases (curve only, no --raw-dir):
  help                 the script starts and prints usage
  curve_sharp          an unambiguous boundary is found
  json_contract        the JSON carries every key the dashboard reads, and the script's own
  deriv_thinned        the derivative arrays are paired and within --export-points
  cutoff_format        cutoff.txt is a bare integer with no trailing newline
  cutoff_positions     the inclusive cutoff keeps at least the centroid's cells
  curve_deterministic  two runs write byte-identical outputs
  manual_cutoff        -m bypasses the search and still exports the curve
  curve_two_knee       the real boundary is taken below a deeper decoy bend that is on offer
  curve_shoulder       the cells are kept, not the ambient shoulder, at 0.5x-5x expected
  curve_no_knee        a library with no boundary warns and is less well localised
  tiny_fallback        a five-barcode curve degrades to the documented fallback

Cases (with --raw-dir, the ambient test):
  fixture_fair                 the curve alone fails the two_tier, shoulder and gradual traps
  shallow_fair                 the curve alone stops at a bend among cells above the floor
  ambient_two_tier             the cell/empty boundary, not a subpopulation bend, with no, 0.1x or 10x expected
  ambient_shoulder             the cells, not the debris shoulder, with no, 0.1x or 10x expected
  ambient_gradual              the gated scan places a boundary with no bend, with no, 0.1x or 10x expected
  ambient_split_pool           with no empty plateau, the cliff is taken and shallow cells kept
  ambient_shallow_split_pool   the cliff below the floor is found (tier 3)
  ambient_json_contract        the JSON carries the ambient test's keys and the selection tier
  ambient_deterministic        two runs write byte-identical outputs
  alevin_input                 an alevin-fry quant directory gives the same cutoff as STARsolo's raw/
  matrix_mismatch              a matrix that does not match the curve is ignored, with a warning
  tier2_z                      a second-tier bend needs a clear contrast z; the most prominent passing one wins
  scan_gate                    the scan is refused when its empties' level is too high to be empty droplets
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
    record SKIP "python" "numpy and scipy are needed by the cell caller"
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
#
# The verdict is read from the last line, so a warning a library prints on
# stderr cannot mask it.
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
        out="$(printf '%s\n' "$out" | tail -1)"
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

# call NAME INPUT PREFIX [EXTRA ARGS...] -- run the cell caller on a fixture.
#
# INPUT is relative to WORKDIR; the outputs are WORKDIR/PREFIX.json and
# WORKDIR/PREFIX_cutoff.txt. A failing run is recorded rather than propagated:
# under `set -e` a non-zero exit here would abort the whole check, losing every
# case after it, when what is wanted is one FAIL and the remaining cases still
# reported.
call() {
    local name="$1" input="$2" prefix="$3"; shift 3
    local logfile="$BCA_TEST_LOGDIR/cell_calling_$(slugify "$prefix").log"
    if ! run_logged "$logfile" \
            "$PYTHON" "$TOOL" \
            -i "$WORKDIR/$input" -s "$name" \
            -o "$WORKDIR/${prefix}.json" -c "$WORKDIR/${prefix}_cutoff.txt" "$@"; then
        record FAIL "run_$prefix" "the cell caller exited non-zero"
        tail_log "$logfile"
    fi
}

# --------------------------------------------------------------------------
# Case: help
# --------------------------------------------------------------------------

if "$PYTHON" "$TOOL" --help >/dev/null 2>&1; then
    record PASS "help" "secondderiv_cellcalling.py prints usage"
else
    record FAIL "help" "secondderiv_cellcalling.py --help failed"
    finish_check
    exit 1
fi

# --------------------------------------------------------------------------
# Curve-only fixtures
# --------------------------------------------------------------------------

if ! run_logged "$BCA_TEST_LOGDIR/cell_calling_fixture.log" \
        "$PYTHON" "$FIXTURE_GEN" --outdir "$WORKDIR"; then
    record FAIL "fixture" "could not build the fixtures"
    tail_log "$BCA_TEST_LOGDIR/cell_calling_fixture.log"
    finish_check
    exit 1
fi
record PASS "fixture" "four archetypes plus a five-barcode curve"

call sharp sharp_UMIperCellSorted.txt curve_sharp -e 800
call sharp sharp_UMIperCellSorted.txt curve_sharp_again -e 800
call sharp sharp_UMIperCellSorted.txt curve_sharp_centroid -e 800 --cutoff-position centroid
call sharp sharp_UMIperCellSorted.txt manual -m 777
call two_knee two_knee_UMIperCellSorted.txt curve_two_knee -e 1000
call no_knee no_knee_UMIperCellSorted.txt curve_no_knee -e 2000
for factor in 0.5 1 5; do
    call ambient_shoulder ambient_shoulder_UMIperCellSorted.txt "curve_shoulder_$factor" \
        -e "$(awk -v f="$factor" 'BEGIN { printf "%d", 2500 * f }')"
done
call tiny tiny_UMIperCellSorted.txt tiny -e 800

# --------------------------------------------------------------------------
# Case: sharp -- an unambiguous boundary, which the curve rule must find
# --------------------------------------------------------------------------

# The inclusive position keeps the low-count end of the bend, so it is allowed to
# land a little past the planted rank; the tolerance band is wide enough for that.
assert "curve_sharp" \
    "(lambda e: (within_tolerance(e.get('cutoff_rank'), 'sharp')[0], \
                 'rank %s vs planted %s (%.3f decades), rule %s' % (e.get('cutoff_rank'), PLANTED_BOUNDARY['sharp'], \
                 within_tolerance(e.get('cutoff_rank'), 'sharp')[1], e.get('selection_rule'))))(entry(A[0]))" \
    "" "$WORKDIR/curve_sharp.json"

# --------------------------------------------------------------------------
# Case: the output contract the dashboard depends on
# --------------------------------------------------------------------------

assert "json_contract" \
    "(lambda e: (all(k in e for k in LEGACY_KEYS + ['script_version', 'selection_rule', 'selection_tier', 'cutoff_position', \
                                                   'cutoff_rank', 'cutoff_rank_centroid', 'cutoff_rank_inclusive', \
                                                   'n_qualifying', 'candidates', 'ambient_status', 'ambient_tests']) \
                 and e['script_version'] == 6 and e['status'] in ('ok', 'warning') \
                 and len(e['logX']) == len(e['logY']) == len(e['customdata']) \
                 and e['ambient_status'] == 'no_matrix' and e['selection_tier'] == 'curve', \
                 'missing: %s; ambient_status %s, tier %s' % ([k for k in LEGACY_KEYS if k not in e], \
                 e.get('ambient_status'), e.get('selection_tier'))))(entry(A[0]))" \
    "every key the report reads, plus the script's own" "$WORKDIR/curve_sharp.json"

assert "deriv_thinned" \
    "(lambda e: (len(e['derivX']) == len(e['derivY']) and 0 < len(e['derivX']) <= 2000 and len(e['raw_slope_derivY']) == len(e['derivX']), \
                 '%d derivative points' % len(e['derivX'])))(entry(A[0]))" \
    "" "$WORKDIR/curve_sharp.json"

assert "cutoff_format" \
    "(lambda b: (b == str(int(b)).encode() and int(b) == entry(A[1])['threshold_umi'], \
                 'cutoff.txt is %r' % b))(raw(A[0]))" \
    "bare integer, no trailing newline" "$WORKDIR/curve_sharp_cutoff.txt" "$WORKDIR/curve_sharp.json"

assert "cutoff_positions" \
    "(lambda inc, cen: (inc['cutoff_position'] == 'inclusive' and cen['cutoff_position'] == 'centroid' \
                        and inc['threshold_umi'] <= cen['threshold_umi'] and cen['threshold_umi'] == inc['cutoff_umi_centroid'], \
                        'inclusive %s UMIs, centroid %s UMIs' % (inc['threshold_umi'], cen['threshold_umi'])))(entry(A[0]), entry(A[1]))" \
    "" "$WORKDIR/curve_sharp.json" "$WORKDIR/curve_sharp_centroid.json"

assert "curve_deterministic" \
    "(raw(A[0]) == raw(A[1]) and raw(A[2]) == raw(A[3]), 'identical' if raw(A[0]) == raw(A[1]) else 'JSON differs between runs')" \
    "two runs, byte-identical JSON and cutoff" \
    "$WORKDIR/curve_sharp.json" "$WORKDIR/curve_sharp_again.json" \
    "$WORKDIR/curve_sharp_cutoff.txt" "$WORKDIR/curve_sharp_again_cutoff.txt"

# --------------------------------------------------------------------------
# Case: manual cutoff -- the search is skipped, the curve is not
# --------------------------------------------------------------------------

assert "manual_cutoff" \
    "(lambda e, b: (e['threshold_umi'] == 777 and int(b) == 777 and len(e['logX']) > 0 and e['selection_rule'] == 'manual', \
                    'threshold_umi=%s, cutoff.txt=%s, rule=%s' % (e['threshold_umi'], b.decode(), e['selection_rule'])))(entry(A[0]), raw(A[1]))" \
    "" "$WORKDIR/manual.json" "$WORKDIR/manual_cutoff.txt"

# --------------------------------------------------------------------------
# Case: two competing bends -- the selection rule is the whole point
# --------------------------------------------------------------------------

# Both halves matter. Landing on the boundary is the claim; a candidate at the
# decoy is what shows the fixture offered the choice. An earlier fixture put the
# decoy where no rule could reach it, and this case passed without exercising
# the choice at all.
assert "curve_two_knee" \
    "(lambda e: (within_tolerance(e.get('cutoff_rank'), 'two_knee')[0] \
                 and any(abs(np.log10(c['rank'] / DECOY_BEND['two_knee'])) <= TOLERANCE_DECADES for c in e['candidates']), \
                 'rank %s (planted %s), rule %s; candidates at the decoy %s: %s' \
                 % (e.get('cutoff_rank'), PLANTED_BOUNDARY['two_knee'], e.get('selection_rule'), DECOY_BEND['two_knee'], \
                    ['rank %d, prominence %.2f' % (c['rank'], c['prominence']) for c in e['candidates'] \
                     if abs(np.log10(c['rank'] / DECOY_BEND['two_knee'])) <= TOLERANCE_DECADES])))(entry(A[0]))" \
    "" "$WORKDIR/curve_two_knee.json"

# The shoulder's lower edge is the deepest bend on the curve, and an over-set
# expected count must not hand it the whole shoulder.
assert "curve_shoulder" \
    "(lambda es: (all(within_tolerance(e.get('cutoff_rank'), 'ambient_shoulder')[0] for e in es), \
                  'ranks %s at 0.5x/1x/5x expected vs planted %s (shoulder edge %s)' \
                  % ([e.get('cutoff_rank') for e in es], PLANTED_BOUNDARY['ambient_shoulder'], DECOY_BEND['ambient_shoulder'])))([entry(a) for a in A])" \
    "" "$WORKDIR/curve_shoulder_0.5.json" "$WORKDIR/curve_shoulder_1.json" "$WORKDIR/curve_shoulder_5.json"

# --------------------------------------------------------------------------
# Case: no boundary at all -- the answer is that there is no answer
# --------------------------------------------------------------------------

# A curve with no boundary must be less well localised than one with an obvious
# boundary, or the resolution diagnostic means nothing, and the run must say it
# found nothing defensible. It still writes a cutoff -- three downstream modules
# read the file -- but under a warning.
assert "curve_no_knee" \
    "(lambda flat, sharp: (flat['status'] == 'warning' \
                           and (flat.get('resolution_decades') or 0) > (sharp.get('resolution_decades') or 0), \
                           'status %s, rule %s; %s decades indistinguishable on a featureless curve vs %s on a sharp one' \
                           % (flat['status'], flat.get('selection_rule'), flat.get('resolution_decades'), \
                              sharp.get('resolution_decades'))))(entry(A[0]), entry(A[1]))" \
    "" "$WORKDIR/curve_no_knee.json" "$WORKDIR/curve_sharp.json"

# --------------------------------------------------------------------------
# Case: too little data -- the documented fallback, not a crash
# --------------------------------------------------------------------------

assert "tiny_fallback" \
    "(lambda e, b: (e['logX'] == [] and e['status'] == 'warning' and e['threshold_umi'] == 0 and int(b) == 0, \
                    'status=%s, threshold=%s' % (e['status'], e['threshold_umi'])))(entry(A[0]), raw(A[1]))" \
    "" "$WORKDIR/tiny.json" "$WORKDIR/tiny_cutoff.txt"

# --------------------------------------------------------------------------
# Matrix fixtures -- the bend chosen by what the barcodes either side look like
# --------------------------------------------------------------------------

# STARsolo raw directories with planted expression (tests/lib/make_ambient_fixture.py):
# a high-RNA subpopulation above the real cells, as on 2026_Q3 spolac01_lib1; an ambient
# debris shoulder below them; cells that thin out with no bend at the boundary, as on
# spolac01_lib2; and two split-pool (Parse) libraries with no empty plateau. Every
# ambient case runs with no expected count, 0.1x and 10x, since the samplesheet's is
# often a placeholder.
if run_logged "$BCA_TEST_LOGDIR/cell_calling_ambient_fixture.log" \
        "$PYTHON" "$AMBIENT_FIXTURE_GEN" --outdir "$WORKDIR"; then
    record PASS "ambient_fixture" "five archetypes with planted expression"

    for kind in split_pool shallow_split_pool two_tier shoulder gradual; do
        D="${kind}_Solo.out/GeneFull_Ex50pAS"
        exp="$("$PYTHON" -c "import sys; sys.path.insert(0, sys.argv[1]); import make_ambient_fixture as m; print(m.EXPECTED_CELLS[sys.argv[2]])" "$TESTS_DIR/lib" "$kind" | tail -1)"
        call "$kind" "$D/UMIperCellSorted.txt" "ambient_${kind}_none" --raw-dir "$WORKDIR/$D/raw"
        call "$kind" "$D/UMIperCellSorted.txt" "ambient_${kind}_low" -e "$((exp / 10))" --raw-dir "$WORKDIR/$D/raw"
        call "$kind" "$D/UMIperCellSorted.txt" "ambient_${kind}_high" -e "$((exp * 10))" --raw-dir "$WORKDIR/$D/raw"
    done
    SP="split_pool_Solo.out/GeneFull_Ex50pAS"
    SSP="shallow_split_pool_Solo.out/GeneFull_Ex50pAS"
    TT="two_tier_Solo.out/GeneFull_Ex50pAS"
    SH="shoulder_Solo.out/GeneFull_Ex50pAS"
    GR="gradual_Solo.out/GeneFull_Ex50pAS"
    call split_pool "$SP/UMIperCellSorted.txt" ambient_split_pool_none_again --raw-dir "$WORKDIR/$SP/raw"

    # The same curves without --raw-dir: the curve rule each trap is built to defeat
    call two_tier "$TT/UMIperCellSorted.txt" curveonly_two_tier_450 -e 450
    call two_tier "$TT/UMIperCellSorted.txt" curveonly_two_tier_3000 -e 3000
    call shoulder "$SH/UMIperCellSorted.txt" curveonly_shoulder_12500 -e 12500
    call gradual "$GR/UMIperCellSorted.txt" curveonly_gradual_4000 -e 4000
    call shallow_split_pool "$SSP/UMIperCellSorted.txt" curveonly_shallow_split_pool_4000 -e 4000

    # The two_tier curve against the shoulder's matrix: the totals disagree rank for rank
    call two_tier "$TT/UMIperCellSorted.txt" matrix_mismatch -e 3000 --raw-dir "$WORKDIR/$SH/raw"

    # F1 against a fixture's planted cells, the best any cutoff reaches, and a rank
    # within TOLERANCE_DECADES of the planted boundary, all named by the generator
    BUILD="(lambda kind: __import__('make_ambient_fixture').build(kind))"
    F1="(lambda cut, d: (lambda pr: 2 * pr[0] * pr[1] / (pr[0] + pr[1]) if sum(pr) else 0.0) \
         (__import__('make_ambient_fixture').precision_recall(d[1], d[2], cut)))"
    BEST_F1="(lambda d: max($F1(c, d) for c in np.unique(d[1]) if c >= 20))"
    NEAR="(lambda e, kind: (lambda b: abs(np.log10(e['cutoff_rank'] / b)) <= TOLERANCE_DECADES)(__import__('make_ambient_fixture').PLANTED_BOUNDARY[kind]))"

    # ----------------------------------------------------------------------
    # Fairness: the curve alone fails each trap, so the ambient cases test something
    # ----------------------------------------------------------------------

    assert "fixture_fair" \
        "(lambda tt, sh, gr, d: (not $NEAR(tt, 'two_tier') and not $NEAR(sh, 'shoulder') and $F1(gr['threshold_umi'], d) < 0.6, \
                                 'curve alone: rank %s on two_tier at 0.15x and %s on shoulder at 5x; F1 %.2f on gradual at 1x' \
                                 % (tt['cutoff_rank'], sh['cutoff_rank'], $F1(gr['threshold_umi'], d)))) \
         (entry(A[0]), entry(A[1]), entry(A[2]), $BUILD('gradual'))" \
        "" "$WORKDIR/curveonly_two_tier_450.json" "$WORKDIR/curveonly_shoulder_12500.json" "$WORKDIR/curveonly_gradual_4000.json"

    # The cells straddle the floor: deep types above it, a 45-UMI type below, so the cliff
    # lies under 100 UMIs and a bend among the cells sits above it, where the curve stops.
    assert "shallow_fair" \
        "(lambda e, d: (lambda best: ($F1(e['threshold_umi'], d) < 0.9 * best, \
                                      'curve alone cuts at %s UMIs (%s), F1 %.3f vs best possible %.3f' \
                                      % (e['threshold_umi'], e['selection_rule'], $F1(e['threshold_umi'], d), best)))($BEST_F1(d))) \
         (entry(A[0]), $BUILD('shallow_split_pool'))" \
        "" "$WORKDIR/curveonly_shallow_split_pool_4000.json"

    # split_pool has no fairness case: the curve alone lands on its cliff too. What it
    # defeats is an ungated ambient scan, and scan_gate below covers that rule directly.

    # ----------------------------------------------------------------------
    # The ambient test places the boundary, whatever the expected count
    # ----------------------------------------------------------------------

    for kind in two_tier shoulder; do
        assert "ambient_$kind" \
            "(lambda es: (all($NEAR(e, '$kind') and e['selection_rule'].startswith('ambient') for e in es), \
                          'ranks %s (%s) with no/0.1x/10x expected vs planted %s, decoy bend %s' \
                          % ([e['cutoff_rank'] for e in es], [e['selection_rule'] for e in es], \
                             __import__('make_ambient_fixture').PLANTED_BOUNDARY['$kind'], \
                             __import__('make_ambient_fixture').DECOY_BEND['$kind'])))([entry(a) for a in A])" \
            "" "$WORKDIR/ambient_${kind}_none.json" "$WORKDIR/ambient_${kind}_low.json" "$WORKDIR/ambient_${kind}_high.json"
    done

    # Cells and empties overlap in depth here, so these are scored against the planted
    # cells: nearly as well as any cutoff could
    for kind in gradual split_pool shallow_split_pool; do
        assert "ambient_$kind" \
            "(lambda es, d: (lambda best: (all(e['selection_rule'].startswith('ambient') and $F1(e['threshold_umi'], d) >= 0.95 * best for e in es), \
                                           'cutoffs %s UMIs (%s), F1 %s with no/0.1x/10x expected vs best possible %.3f' \
                                           % ([e['threshold_umi'] for e in es], [e['selection_rule'] for e in es], \
                                              ['%.3f' % $F1(e['threshold_umi'], d) for e in es], best)))($BEST_F1(d))) \
             ([entry(a) for a in A], $BUILD('$kind'))" \
            "" "$WORKDIR/ambient_${kind}_none.json" "$WORKDIR/ambient_${kind}_low.json" "$WORKDIR/ambient_${kind}_high.json"
    done

    assert "ambient_json_contract" \
        "(lambda e: (all(k in e for k in LEGACY_KEYS + ['script_version', 'selection_rule', 'selection_tier', 'cutoff_rank', \
                                                       'n_qualifying', 'candidates', 'ambient_status', 'ambient_model', 'ambient_tests']) \
                     and e['script_version'] == 6 and e['ambient_status'] == 'ok' and e['selection_tier'] == 2 \
                     and all(0 <= t['f_above'] <= 1 and 0 <= t['f_below'] <= 1 and 'z' in t and t['tier'] in (1, 2) \
                             for t in e['ambient_tests']), \
                     'status %s, tier %s, %d ambient tests' % (e.get('ambient_status'), e.get('selection_tier'), \
                                                                len(e.get('ambient_tests', [])))))(entry(A[0]))" \
        "the dashboard's keys plus the ambient test's and the tier" "$WORKDIR/ambient_split_pool_none.json"

    assert "ambient_deterministic" \
        "(raw(A[0]) == raw(A[1]) and raw(A[2]) == raw(A[3]), 'identical' if raw(A[0]) == raw(A[1]) else 'JSON differs between runs')" \
        "two runs, byte-identical JSON and cutoff" \
        "$WORKDIR/ambient_split_pool_none.json" "$WORKDIR/ambient_split_pool_none_again.json" \
        "$WORKDIR/ambient_split_pool_none_cutoff.txt" "$WORKDIR/ambient_split_pool_none_again_cutoff.txt"

    # A raw directory the curve was not counted from must not steer the choice: the
    # ambient test is skipped with a warning, and the cutoff is the curve rule's.
    assert "matrix_mismatch" \
        "(lambda e, b, ref: (e['ambient_status'] == 'no_matrix' and 'does not match' in e['message'] and b == ref, \
                             'ambient_status %s, cutoff %s vs %s from the curve alone' \
                             % (e['ambient_status'], b.decode(), ref.decode())))(entry(A[0]), raw(A[1]), raw(A[2]))" \
        "" "$WORKDIR/matrix_mismatch.json" "$WORKDIR/matrix_mismatch_cutoff.txt" "$WORKDIR/curveonly_two_tier_3000_cutoff.txt"

    # The same library as an alevin-fry quant directory: cells x USA columns (the bare gene
    # IDs, then <gene>-U, then <gene>-A), the counts split over the S and U blocks. The
    # curve is derived as SECONDDERIV_CELLCALLING derives it, and the caller must read the
    # same counts from the directory and land on the same cutoff.
    AF="$WORKDIR/split_pool_alevin"
    if run_logged "$BCA_TEST_LOGDIR/cell_calling_alevin_fixture.log" "$PYTHON" -c "
import os, sys
import numpy as np, scipy.io, scipy.sparse as sp
sys.path.insert(0, sys.argv[1])
from make_ambient_fixture import build
counts, _totals, _cells = build('split_pool')
u = counts // 3
usa = np.hstack([counts - u, u, np.zeros_like(counts)])
out = sys.argv[2]
os.makedirs(out, exist_ok=True)
scipy.io.mmwrite(os.path.join(out, 'quants_mat.mtx'), sp.coo_matrix(usa), field='integer')
genes = ['g%d' % g for g in range(counts.shape[1])]
with open(os.path.join(out, 'quants_mat_rows.txt'), 'w') as fh:
    fh.write('\n'.join('BC%06d' % i for i in range(counts.shape[0])) + '\n')
with open(os.path.join(out, 'quants_mat_cols.txt'), 'w') as fh:
    fh.write('\n'.join(genes + [g + '-U' for g in genes] + [g + '-A' for g in genes]) + '\n')
" "$TESTS_DIR/lib" "$AF" \
        && run_logged "$BCA_TEST_LOGDIR/cell_calling_alevin_curve.log" "$PYTHON" "$ALEVIN_TOOL" \
            umis -d "$AF" -o "$WORKDIR/split_pool_alevin_UMIperCellSorted.txt" --counts SUA; then
        call split_pool split_pool_alevin_UMIperCellSorted.txt split_pool_alevin --raw-dir "$AF" --counts SUA
        assert "alevin_input" \
            "(lambda a, s, e: (raw(a) == raw(s) and e['ambient_status'] == 'ok', \
                              'cutoff %s UMIs from the alevin-fry directory, %s from STARsolo raw/; ambient test %s' \
                              % (raw(a).decode(), raw(s).decode(), e['ambient_status'])))(A[0], A[1], entry(A[2]))" \
            "" "$WORKDIR/split_pool_alevin_cutoff.txt" "$WORKDIR/ambient_split_pool_none_cutoff.txt" \
            "$WORKDIR/split_pool_alevin.json"
    else
        record FAIL "alevin_input" "could not build the alevin-fry fixture or its curve"
        tail_log "$BCA_TEST_LOGDIR/cell_calling_alevin_fixture.log"
    fi
else
    record FAIL "ambient_fixture" "could not build the matrix fixtures"
    tail_log "$BCA_TEST_LOGDIR/cell_calling_ambient_fixture.log"
fi

# --------------------------------------------------------------------------
# The selection rules themselves, with the ambient test stubbed
# --------------------------------------------------------------------------

# A second-tier bend needs a clear z, the most prominent passing one wins, a
# qualifying boundary comes first, and the scan is refused when its "empties" are
# too often non-ambient.
if run_logged "$BCA_TEST_LOGDIR/cell_calling_rules.log" "$PYTHON" -c "
import importlib.util, json, sys
import numpy as np
spec = importlib.util.spec_from_file_location('cellcalling', sys.argv[1])
cc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cc)
bands = {}
def band(model, matrix, rank, seed):
    fa, na, fb, nb = bands[rank]
    return {'f_above': fa, 'n_above': na, 'f_below': fb, 'n_below': nb, 'contrast': fa - fb,
            'z': cc.contrast_z(fa, na, fb, nb)}
cc.band_contrast = band
def cand(rank, prominence, local, qualifies=False):
    return {'rank': rank, 'index': rank - 1, 'prominence': prominence, 'local_sigmas': local,
            'qualifies': qualifies, 'global_sigmas': 1.0, 'recovery_sigmas': 0.0}
out = {}
bands = {1000: (0.60, 40, 0.50, 40)}
out['small_bands'] = cc.select_by_ambient([cand(1000, 2.0, 20.0)], None, None, None, 0)[1]
bands = {1000: (0.9, 2000, 0.3, 2000), 3000: (0.8, 2000, 0.2, 2000)}
c, how = cc.select_by_ambient([cand(1000, 3.0, 20.0), cand(3000, 2.0, 20.0)], None, None, None, 0)
out['most_prominent'] = [how, c['rank']]
c, how = cc.select_by_ambient([cand(1000, 3.0, 20.0), cand(3000, 2.0, 20.0, qualifies=True)], None, None, None, 0)
out['tier1_first'] = [how, c['rank']]
out['insignificant'] = cc.select_by_ambient([cand(1000, 3.0, 4.0)], None, None, None, 0)[1]
# Scan: ranks below 1000 are cells, the rest non-ambient at a given level
def stub_pvalues(level):
    return lambda model, X, rows, rng: np.where((rows < 1000) | ((rows % 100) < level * 100), 0.0, 1.0)
matrix = {'X': None, 'totals': np.arange(5000, 0, -1).astype(float)}
for name, level in (('scan_shallow_cells', 0.3), ('scan_empties', 0.02)):
    cc.ambient_pvalues = stub_pvalues(level)
    s = cc.ambient_scan(None, matrix, 5000, 0)
    out[name] = [s['gate'], s['rank'], round(s['empty_level'], 3)]
json.dump(out, open(sys.argv[2], 'w'))
" "$TOOL" "$WORKDIR/rules.json"; then
    assert "tier2_z" \
        "(lambda r: (r['small_bands'] == 'no_ambient_boundary' and r['most_prominent'] == ['ambient_tier2', 1000] \
                     and r['tier1_first'] == ['ambient_contrast', 3000] and r['insignificant'] == 'no_ambient_boundary', \
                     str(r)))(json.load(open(A[0])))" \
        "" "$WORKDIR/rules.json"
    assert "scan_gate" \
        "(lambda r: (r['scan_shallow_cells'][0] == 'empties_not_reached' and r['scan_shallow_cells'][1] is None \
                     and r['scan_empties'][0] == 'passed' and abs(np.log10(r['scan_empties'][1] / 1000)) <= 0.1, \
                     'shallow cells %s, empties %s' % (r['scan_shallow_cells'], r['scan_empties'])))(json.load(open(A[0])))" \
        "" "$WORKDIR/rules.json"
else
    record FAIL "rules" "the stubbed selection-rule script failed"
    tail_log "$BCA_TEST_LOGDIR/cell_calling_rules.log"
fi

finish_check
