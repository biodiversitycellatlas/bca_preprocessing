#!/usr/bin/env bash
# description: Parse synthetic Nextflow traces and validate the tuning recommendations
#
# bin/resource_efficiency.py turns execution traces into scheduled-resource
# recommendations, so a silent parsing error there does not raise an exception --
# it produces a plausible-looking number that is wrong, and the next run is sized
# from it. Nothing on a developer machine has a real trace to check against, so
# this builds synthetic ones with a known scaling exponent (tests/lib/
# make_trace_fixture.py) plus a hand-written file covering the awkward formats
# Nextflow emits (tests/fixtures/execution_trace_formats.txt), and asserts on the
# numbers that come back out.
#
# No sequencing data and no cluster are needed.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

TOOL="$PROJECT_ROOT/bin/resource_efficiency.py"
FIXTURE_GEN="$TESTS_DIR/lib/make_trace_fixture.py"
FORMAT_FIXTURE="$TESTS_DIR/fixtures/execution_trace_formats.txt"

PYTHON="${BCA_PYTHON:-}"
EXPONENT="0.85"
TOLERANCE="0.15"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh resource_efficiency [-- OPTIONS]
       tests/checks/resource_efficiency.sh [OPTIONS]

Validate bin/resource_efficiency.py against synthetic execution traces.

Options:
  --python PATH       Python interpreter to use (default: \$BCA_PYTHON, else
                      python3, else python).
  --exponent FLOAT    Scaling exponent to embed in the generated runs.
                      Default: ${EXPONENT}.
  --tolerance FLOAT   Allowed absolute error when recovering it.
                      Default: ${TOLERANCE}.
  --keep              Keep the generated fixtures for inspection. They are
                      written under \$BCA_TEST_LOGDIR either way; this only
                      stops them being deleted at the end.
  -h, --help          Show this message.

Cases:
  help                    the tool starts and prints usage
  trace_formats           byte/duration/percent parsing, CACHED and FAILED
                          handling, aliases, unlabelled and unknown processes
  label_coverage          every process in modules/ maps to a tier
  base_config             conf/base.config tiers parse to their declared values
  scaling_recovery        the embedded exponent comes back within tolerance
  recommendation_sanity   no recommendation is below what was observed
  single_run_fallback     one run yields flat recommendations, never a fit
  coverage_guard          a partly observed tier is never lowered
  emitted_config          the generated config is valid Nextflow (needs nextflow)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --python)    PYTHON="${2:?--python needs a value}"; shift 2 ;;
        --exponent)  EXPONENT="${2:?--exponent needs a value}"; shift 2 ;;
        --tolerance) TOLERANCE="${2:?--tolerance needs a value}"; shift 2 ;;
        --keep)      KEEP=1; shift ;;
        -h|--help)   usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

# --------------------------------------------------------------------------
# Pre-flight
# --------------------------------------------------------------------------

log_header "resource_efficiency"

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

WORKDIR="$BCA_TEST_LOGDIR/resource_efficiency"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# jq is not assumed; the JSON assertions are made in python, which is already required.
assert_json() {
    # assert_json NAME JSON_FILE PYTHON_EXPR DESCRIPTION
    local name="$1" json="$2" expr="$3" detail="${4:-}"
    local out
    if out="$("$PYTHON" -c "
import json, sys
d = json.load(open(sys.argv[1]))
procs = {p['process']: p for p in d['processes']}
ok, msg = ($expr)
print(('PASS' if ok else 'FAIL') + '\t' + str(msg))
" "$json" 2>&1)"; then
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

if run_logged "$BCA_TEST_LOGDIR/resource_efficiency_help.log" "$PYTHON" "$TOOL" --help; then
    record PASS "help" "--help exits 0"
else
    record FAIL "help" "--help failed, see resource_efficiency_help.log"
    finish_check
    exit 1
fi

# --------------------------------------------------------------------------
# Case: trace_formats
#
# The hand-written fixture is the only ground truth available for the parsing
# layer, so it is checked before anything statistical is trusted.
# --------------------------------------------------------------------------

FMT_DIR="$WORKDIR/formats"
mkdir -p "$FMT_DIR/pipeline_info"
cp "$FORMAT_FIXTURE" "$FMT_DIR/pipeline_info/execution_trace_2026-03-01_09-00-00.txt"

FMT_JSON="$WORKDIR/formats.json"
if run_logged "$BCA_TEST_LOGDIR/resource_efficiency_formats.log" \
        "$PYTHON" "$TOOL" --results "$FMT_DIR" --no-plots --min-tasks 1 \
        --output "$WORKDIR/formats_report" --json "$FMT_JSON"; then
    record PASS "trace_formats.run" "parsed the format fixture"
else
    record FAIL "trace_formats.run" "see resource_efficiency_formats.log"
    finish_check
    exit 1
fi

# 15 data rows: 11 COMPLETED, 2 FAILED, 1 CACHED, 1 ABORTED.
assert_json "trace_formats.counts" "$FMT_JSON" \
    "(d['n_completed_tasks'] == 11 and d['n_failed_tasks'] == 2 and d['n_cached_rows'] == 1,
      f\"completed={d['n_completed_tasks']} failed={d['n_failed_tasks']} cached={d['n_cached_rows']}\")" \
    "11 completed, 2 failed, 1 cached skipped"

# 1.2 GB must be read base-1024. Reading it as 1.2e9 understates memory by 7%
# on every single row, which is exactly the kind of error that looks plausible.
assert_json "trace_formats.base1024" "$FMT_JSON" \
    "(abs(procs['FASTQC']['max_peak_rss'] - 2.4 * 1024**3) < 1,
      f\"FASTQC peak_rss={procs['FASTQC']['max_peak_rss']}, expected {2.4 * 1024**3}\")" \
    "2.4 GB parsed as 2.4 GiB"

# 1d 2h == 93600 s. A compound duration with a day component is easy to drop.
assert_json "trace_formats.duration" "$FMT_JSON" \
    "(abs(procs['STARSOLO_INDEX']['p95_realtime'] - 93600) < 1,
      f\"STARSOLO_INDEX realtime={procs['STARSOLO_INDEX']['p95_realtime']}, expected 93600\")" \
    "'1d 2h' parsed as 93600 s"

# The exit-137 row must be recognised as a kill and force the recommendation up;
# the exit-1 row is an ordinary tool error and must not.
assert_json "trace_formats.killed" "$FMT_JSON" \
    "(sorted((f['process'], f['killed']) for f in d['failures'])
        == [('SCRUBLET', False), ('STARSOLO_ALIGN', True)],
      str(sorted((f['process'], f['killed']) for f in d['failures'])))" \
    "exit 137 is a kill, exit 1 is not"

# Aliased inclusions are recorded under the alias but declare no label of their own.
assert_json "trace_formats.aliases" "$FMT_JSON" \
    "(procs['DOUBLET_FILTER_RAW']['label'] == 'process_low'
      and procs['MERGE_REF_GTF_GENEEXT']['label'] == '__default__',
      f\"DOUBLET_FILTER_RAW={procs['DOUBLET_FILTER_RAW']['label']} \"
      f\"MERGE_REF_GTF_GENEEXT={procs['MERGE_REF_GTF_GENEEXT']['label']}\")" \
    "aliases resolve to their module's tier"

# A process no module declares must be surfaced, never silently folded into a tier.
assert_json "trace_formats.unmapped" "$FMT_JSON" \
    "(d['unmapped_processes'] == ['NOT_A_REAL_PROCESS'], str(d['unmapped_processes']))" \
    "unknown process reported, not swallowed"

# A row with '-' in every usage column must degrade that metric only.
assert_json "trace_formats.missing_values" "$FMT_JSON" \
    "(procs['MTX_TO_10X']['max_peak_rss'] is None
      and procs['MTX_TO_10X']['p95_realtime'] is not None,
      f\"peak_rss={procs['MTX_TO_10X']['max_peak_rss']} realtime={procs['MTX_TO_10X']['p95_realtime']}\")" \
    "'-' yields None without losing the row"

# --------------------------------------------------------------------------
# Case: label_coverage and base_config
#
# Guards against a new module whose label the mapper cannot see, and against the
# base.config parser silently returning None for a tier.
# --------------------------------------------------------------------------

COVERAGE_OUT="$WORKDIR/coverage.txt"
if "$PYTHON" - "$PROJECT_ROOT" >"$COVERAGE_OUT" 2>&1 <<'PYEOF'
import importlib.util
import os
import sys

root = sys.argv[1]
spec = importlib.util.spec_from_file_location(
    "re_mod", os.path.join(root, "bin", "resource_efficiency.py"))
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

labels = mod.scan_module_labels(root)
if not labels:
    print("FAIL\tno processes found under modules/")
    raise SystemExit(0)

known = set(mod.RESOURCE_LABELS) | {mod.DEFAULT_TIER}
bad = sorted(name for name, info in labels.items() if info.tier not in known)
if bad:
    print("FAIL\tprocesses with an unrecognised tier: " + ", ".join(bad))
    raise SystemExit(0)

# Aliases the workflows declare must resolve back to a real module process.
aliases = ["DOUBLET_FILTER_RAW", "DOUBLET_FILTER_CELL_CALLED",
           "MERGE_REF_GTF_GENEEXT", "MERGE_REF_FASTA_GENEEXT"]
unresolved = [a for a in aliases if mod.resolve_label(a, labels)[1] is None]
if unresolved:
    print("FAIL\taliases did not resolve: " + ", ".join(unresolved))
    raise SystemExit(0)

print(f"PASS\t{len(labels)} processes, all mapped; {len(aliases)} aliases resolved")

# Every tier that declares resources in base.config must parse to a real value.
tiers = mod.parse_base_config(os.path.join(root, "conf", "base.config"))
missing = []
for label in ("process_single", "process_low", "process_medium",
              "process_high", "process_high_memory"):
    spec_ = tiers.get(label)
    if spec_ is None or spec_.cpus is None or spec_.memory is None or spec_.time is None:
        missing.append(label)
if missing:
    print("FAIL\ttiers failed to parse: " + ", ".join(missing))
else:
    summary = ", ".join(
        f"{label}={int(tiers[label].memory / 1024 ** 3)}GB"
        for label in ("process_low", "process_medium", "process_high_memory"))
    print(f"PASS\tall tiers parsed ({summary})")
PYEOF
then
    n=0
    while IFS=$'\t' read -r status detail; do
        n=$((n + 1))
        case "$n" in
            1) name="label_coverage" ;;
            2) name="base_config" ;;
            *) name="label_extra_$n" ;;
        esac
        record "$status" "$name" "$detail"
    done <"$COVERAGE_OUT"
    [[ "$n" -ge 2 ]] || record FAIL "base_config" "no result produced"
else
    record FAIL "label_coverage" "$(tail -3 "$COVERAGE_OUT" | tr '\n' ' ')"
fi

# --------------------------------------------------------------------------
# Case: scaling_recovery, recommendation_sanity, coverage_guard
# --------------------------------------------------------------------------

MULTI_DIR="$WORKDIR/multi"
MULTI_JSON="$WORKDIR/multi.json"
if run_logged "$BCA_TEST_LOGDIR/resource_efficiency_fixture.log" \
        "$PYTHON" "$FIXTURE_GEN" --outdir "$MULTI_DIR" --runs 5 --samples 4 \
        --exponent "$EXPONENT" --seed 1; then
    record PASS "fixture.generate" "5 runs generated with exponent $EXPONENT"
else
    record FAIL "fixture.generate" "see resource_efficiency_fixture.log"
    finish_check
    exit 1
fi

if run_logged "$BCA_TEST_LOGDIR/resource_efficiency_multi.log" \
        "$PYTHON" "$TOOL" --results "$MULTI_DIR" --no-plots \
        --output "$WORKDIR/multi_report" --json "$MULTI_JSON"; then
    record PASS "multi_run.run" "analysed 5 runs"
else
    record FAIL "multi_run.run" "see resource_efficiency_multi.log"
    finish_check
    exit 1
fi

# The fixture's memory model carries a constant term, which damps the recovered
# exponent slightly -- hence a tolerance rather than an equality.
assert_json "scaling_recovery" "$MULTI_JSON" \
    "(lambda fit: (fit is not None and fit['trusted']
                   and abs(fit['exponent'] - $EXPONENT) <= $TOLERANCE
                   and fit['r2'] > 0.9,
                   f\"b={fit and round(fit['exponent'], 3)} r2={fit and round(fit['r2'], 3)}\")
     )(procs['STARSOLO_ALIGN']['fits']['memory'])" \
    ""

# Recommending less than a task actually used guarantees an OOM next run.
assert_json "recommendation_sanity.memory" "$MULTI_JSON" \
    "(all(p['recommended_memory'] >= p['max_peak_rss']
          for p in d['processes']
          if p['recommended_memory'] and p['max_peak_rss']),
      'every recommendation covers the observed maximum')" \
    ""

# CPUs are the one resource never recommended upward: over-requesting only wastes
# allocation, so the tool may shrink but must never inflate it.
assert_json "recommendation_sanity.cpus" "$MULTI_JSON" \
    "(all(p['recommended_cpus'] <= p['current_cpus']
          for p in d['processes']
          if p['recommended_cpus'] and p['current_cpus']),
      'no cpu recommendation exceeds the current request')" \
    ""

# The killed task must drive its process above the request that died.
assert_json "recommendation_sanity.killed" "$MULTI_JSON" \
    "(procs['STARSOLO_ALIGN']['forced_by_failure']
      and procs['STARSOLO_ALIGN']['recommended_memory'] > 128 * 1024**3,
      f\"forced={procs['STARSOLO_ALIGN']['forced_by_failure']} \"
      f\"rec={procs['STARSOLO_ALIGN']['recommended_memory']}\")" \
    "kill raises the recommendation above the failed request"

# A tier is only partly exercised by any real run, so it must never be shrunk on
# that evidence -- the members that did not run would be under-provisioned.
MULTI_CONFIG="$(ls -1 "$WORKDIR/multi_report"/resources_tuned_*.config 2>/dev/null | head -1)"
if [[ -n "$MULTI_CONFIG" ]] && grep -q "were not exercised by these runs" "$MULTI_CONFIG"; then
    if grep -q "held at the conf/base.config value rather than lowered" "$MULTI_CONFIG"; then
        record PASS "coverage_guard" "partly observed tiers are not lowered"
    else
        record FAIL "coverage_guard" "partial coverage noted but the value was still lowered"
    fi
else
    record FAIL "coverage_guard" "no partial-coverage annotation in the emitted config"
fi

# --------------------------------------------------------------------------
# Case: single_run_fallback
#
# One run cannot distinguish scaling from sample spread, so no fit may be
# trusted no matter how clean the correlation looks.
# --------------------------------------------------------------------------

SINGLE_DIR="$WORKDIR/single"
SINGLE_JSON="$WORKDIR/single.json"
"$PYTHON" "$FIXTURE_GEN" --outdir "$SINGLE_DIR" --runs 1 --samples 6 \
    --exponent "$EXPONENT" --seed 2 --no-kill \
    >"$BCA_TEST_LOGDIR/resource_efficiency_single_fixture.log" 2>&1

if run_logged "$BCA_TEST_LOGDIR/resource_efficiency_single.log" \
        "$PYTHON" "$TOOL" --results "$SINGLE_DIR" --no-plots \
        --output "$WORKDIR/single_report" --json "$SINGLE_JSON"; then
    assert_json "single_run_fallback" "$SINGLE_JSON" \
        "(all(not (p['fits'].get('memory') or {}).get('trusted', False)
              for p in d['processes'])
          and all('fit@' not in p['basis'] for p in d['processes']),
          'no fit trusted from a single run')" \
        ""
else
    record FAIL "single_run_fallback" "see resource_efficiency_single.log"
fi

# --------------------------------------------------------------------------
# Case: emitted_config
#
# A config the tool cannot itself validate is a config that breaks the next run.
# --------------------------------------------------------------------------

if ! have_cmd nextflow; then
    record SKIP "emitted_config" "nextflow not installed"
elif [[ -z "$MULTI_CONFIG" ]]; then
    record FAIL "emitted_config" "no config was emitted"
else
    if run_logged "$BCA_TEST_LOGDIR/resource_efficiency_config.log" \
            nextflow config -c "$MULTI_CONFIG" "$PROJECT_ROOT"; then
        record PASS "emitted_config" "nextflow config parses the generated file"
    else
        record FAIL "emitted_config" "see resource_efficiency_config.log"
    fi
fi

finish_check
