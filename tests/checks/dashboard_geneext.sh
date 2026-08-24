#!/usr/bin/env bash
# description: Check that GeneExt statistics reach the dashboard's gene-extension tab
#
# The dashboard summarises GeneExt by re-reading the JSON payload GeneExt embeds
# in its own HTML report, falling back to the plain-text log when that report is
# missing. Neither failure mode is loud: a payload that cannot be located yields
# an empty object, which simply hides the tab, and a fallback that matches
# nothing yields a tab full of blanks. This builds a synthetic GeneExt output
# directory (tests/lib/make_geneext_fixture.py) with known numbers and asserts
# they come back out of the generated dashboard.
#
# It also pins the two reductions the dashboard makes: GeneExt's per-gene
# extension table and orphan-peak BED are deliberately not carried over, since
# they dominate the report's size and are one click away in the report itself.
#
# No sequencing data, no GeneExt installation and no cluster are needed.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

TOOL="$PROJECT_ROOT/bin/generate_dashboard.py"
TEMPLATE="$PROJECT_ROOT/bin/dashboard_report.html"
FIXTURE_GEN="$TESTS_DIR/lib/make_geneext_fixture.py"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh dashboard_geneext [-- OPTIONS]
       tests/checks/dashboard_geneext.sh [OPTIONS]

Validate the gene-extension payload bin/generate_dashboard.py embeds.

Options:
  --python PATH   Python interpreter to use (default: \$BCA_PYTHON, else
                  python3, else python).
  --keep          Keep the generated fixtures and dashboards for inspection.
  -h, --help      Show this message.

Cases:
  help              the tool starts and prints usage
  report_payload    GeneExt's numbers survive into the embedded payload
  payload_trimmed   the per-gene table and orphan BED are left behind
  log_fallback      headline numbers are recovered from the log alone
  absent            a run without GeneExt embeds an empty payload
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

log_header "dashboard_geneext"

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

WORKDIR="$BCA_TEST_LOGDIR/dashboard_geneext"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

printf 'analytical_id,base_id,source\nsampleA_geneext_starsolo,sampleA,starsolo\n' \
    >"$WORKDIR/manifest.csv"
printf 'sample,expected_cells\nsampleA,5000\n' >"$WORKDIR/samplesheet.csv"

# assert_payload NAME DASHBOARD_HTML PYTHON_EXPR [DETAIL]
#
# The payload is a JSON block inside the generated HTML, so it is pulled back out
# the same way a browser would find it: by its script element id. PYTHON_EXPR is
# evaluated with `d` bound to the decoded object and must yield (ok, message).
assert_payload() {
    local name="$1" html="$2" expr="$3" detail="${4:-}"
    local out
    if out="$("$PYTHON" -c "
import json, re, sys
html = open(sys.argv[1], encoding='utf-8').read()
m = re.search(r'<script id=\"geneext-data\" type=\"application/json\">(.*?)</script>', html, re.S)
if not m:
    print('FAIL\tno geneext-data block in the generated dashboard')
    raise SystemExit(0)
d = json.loads(m.group(1))
ok, msg = ($expr)
print(('PASS' if ok else 'FAIL') + '\t' + str(msg))
" "$html" 2>&1)"; then
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

# generate OUTPUT_HTML LOGNAME EXTRA_ARGS...
generate() {
    local output="$1" logname="$2"; shift 2
    run_logged "$BCA_TEST_LOGDIR/dashboard_geneext_${logname}.log" \
        "$PYTHON" "$TOOL" \
        --template "$TEMPLATE" \
        --output "$output" \
        --samplesheet "$WORKDIR/samplesheet.csv" \
        --analytical_samples "$WORKDIR/manifest.csv" \
        "$@"
}

# --------------------------------------------------------------------------
# Case: help
# --------------------------------------------------------------------------

if run_logged "$BCA_TEST_LOGDIR/dashboard_geneext_help.log" "$PYTHON" "$TOOL" --help; then
    record PASS "help" "--help exits 0"
else
    record FAIL "help" "--help failed, see dashboard_geneext_help.log"
    finish_check
    exit 1
fi

# --------------------------------------------------------------------------
# Fixture
# --------------------------------------------------------------------------

GENE_EXT_DIR="$WORKDIR/results/gene_ext"
if ! run_logged "$BCA_TEST_LOGDIR/dashboard_geneext_fixture.log" \
        "$PYTHON" "$FIXTURE_GEN" "$GENE_EXT_DIR"; then
    record FAIL "fixture" "see dashboard_geneext_fixture.log"
    finish_check
    exit 1
fi
REPORT="$GENE_EXT_DIR/geneext.gtf.Report.html"
LOG="$GENE_EXT_DIR/geneext.gtf.GeneExt.log"
record PASS "fixture" "synthetic GeneExt report and log written"

# --------------------------------------------------------------------------
# Case: report_payload
#
# The report is the preferred source, so its numbers -- not the log's -- must be
# the ones that come through. The gene count is what tells the two apart.
# --------------------------------------------------------------------------

FULL="$WORKDIR/dashboard_full.html"
if generate "$FULL" "full" --geneext_report "$REPORT" --geneext_log "$LOG"; then
    record PASS "report_payload.run" "dashboard generated with GeneExt inputs"
else
    record FAIL "report_payload.run" "see dashboard_geneext_full.log"
    finish_check
    exit 1
fi

assert_payload "report_payload.source" "$FULL" \
    "(d.get('source') == 'report', f\"source={d.get('source')}\")" \
    "read from the HTML report, not the log"

assert_payload "report_payload.summary" "$FULL" \
    "(d['summary']['n_extended'] == 13418 and d['summary']['n_genes'] == 31949
      and d['summary']['median_ext'] == 1383.0,
      f\"n_extended={d['summary'].get('n_extended')} n_genes={d['summary'].get('n_genes')} \"
      f\"median_ext={d['summary'].get('median_ext')}\")" \
    "extension counts carried over"

# Both histograms are redrawn by the dashboard, so a missing one is a blank plot.
assert_payload "report_payload.histograms" "$FULL" \
    "(len(d['ext_hist']['counts']) == 100 and len(d['cov_hist']['counts_genic']) == 50
      and d['cov_hist']['log10_threshold'] == 0.301,
      f\"ext_bins={len(d.get('ext_hist', {}).get('counts', []))} \"
      f\"cov_bins={len(d.get('cov_hist', {}).get('counts_genic', []))}\")" \
    "extension and coverage distributions present"

assert_payload "report_payload.peak_flow" "$FULL" \
    "(d['peak_flow']['initial_called'] == 244298
      and d['peak_flow']['passed_filtering'] == 128666
      and d['peak_flow']['assigned_to_genes'] == 13418,
      str(d.get('peak_flow')))" \
    "peak-filtering flow carried over"

# --maxdist drives how far any gene can be extended, so the tab states it; it
# lives in fix_info, one level below where the rest of the summary is read from.
assert_payload "report_payload.extension_param" "$FULL" \
    "(d['extension_param']['mode'] == 'auto'
      and d['extension_param']['effective_value_bp'] == 2982,
      str(d.get('extension_param')))" \
    "--maxdist mode and value carried over"

# Only the genome-fix steps that ran are worth reporting; the skipped ones must
# not appear, or the tab claims fixes that never happened.
assert_payload "report_payload.genome_fix" "$FULL" \
    "(list(d.get('genome_fix', {})) == ['gene_features_added'],
      str(list(d.get('genome_fix', {}))))" \
    "only applied fix steps kept"

# --------------------------------------------------------------------------
# Case: payload_trimmed
#
# GeneExt's report embeds one row per extended gene plus the orphan-peak BED.
# Carrying them into the dashboard would add megabytes per run for data the
# dashboard never draws, so both must be dropped.
# --------------------------------------------------------------------------

assert_payload "payload_trimmed.blocks" "$FULL" \
    "('ext_table' not in d and 'orphan_bed' not in d, str(sorted(d)))" \
    "per-gene table and orphan BED left behind"

assert_payload "payload_trimmed.size" "$FULL" \
    "(len(json.dumps(d)) < 100_000, f\"payload is {len(json.dumps(d))} bytes\")" \
    "payload stays well under the report's own size"

# --------------------------------------------------------------------------
# Case: log_fallback
#
# A GeneExt older than the one that writes an HTML report still has to produce a
# usable tab. The log reports one fewer extended gene than the report does, which
# is how this case is told apart from the one above.
# --------------------------------------------------------------------------

LOGONLY="$WORKDIR/dashboard_logonly.html"
if generate "$LOGONLY" "logonly" --geneext_log "$LOG"; then
    record PASS "log_fallback.run" "dashboard generated from the log alone"
else
    record FAIL "log_fallback.run" "see dashboard_geneext_logonly.log"
    finish_check
    exit 1
fi

assert_payload "log_fallback.numbers" "$LOGONLY" \
    "(d.get('source') == 'log' and d['summary']['n_extended'] == 13417
      and d['summary']['n_genes'] == 31949 and d['summary']['median_ext'] == 1383.0
      and d['peak_flow']['passed_filtering'] == 128666,
      f\"source={d.get('source')} summary={d.get('summary')}\")" \
    "headline numbers recovered from the log"

# --------------------------------------------------------------------------
# Case: absent
#
# Most runs do not set perform_geneext. The payload must then be empty, which is
# what keeps the tab hidden -- not a half-populated object.
# --------------------------------------------------------------------------

NONE="$WORKDIR/dashboard_none.html"
if generate "$NONE" "none"; then
    record PASS "absent.run" "dashboard generated without GeneExt inputs"
else
    record FAIL "absent.run" "see dashboard_geneext_none.log"
    finish_check
    exit 1
fi

assert_payload "absent.empty" "$NONE" \
    "(d == {}, f\"payload={d}\")" \
    "empty payload hides the tab"

finish_check
