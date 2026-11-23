#!/usr/bin/env bash
# description: Check the startup checks on vendor pipeline settings (lib/VendorPipelines.groovy)
#
# A vendor pipeline (Cell Ranger, split-pipe, BD Rhapsody) only runs when it is enabled.
# Every failure here used to be silent: an 'external_pipeline_only' run with nothing to
# run pre-processed the reads and reported success, a split-pipe path that does not exist
# skipped split-pipe, and settings for another protocol were ignored. Each case gives
# VendorPipelines.check a params map and compares its errors and warnings with the
# expected ones, with the container engine of the run where it matters.
#
# Needs nextflow, and is skipped without it.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh vendor_pipelines [-- OPTIONS]
       tests/checks/vendor_pipelines.sh [OPTIONS]

Check lib/VendorPipelines.groovy.

Options:
  --keep          Keep the generated project and outputs for inspection.
  -h, --help      Show this message.

Cases:
  external_only   external_pipeline_only stops when no vendor pipeline would run
  paths           a vendor path that does not exist, or split-pipe without its conda env, stops the run
  mismatch        settings of another protocol's vendor pipeline only warn
  post_mapping    vendor settings only warn, even a path that does not exist
  engine          a vendor pipeline under a profile it cannot run under warns
  clean           defaults, a complete setup and an unknown protocol give nothing
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep)    KEEP=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

log_header "vendor_pipelines"

if ! have_cmd nextflow; then
    record SKIP "nextflow" "nextflow not found"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/vendor_pipelines"
rm -rf "$WORKDIR"
mkdir -p "$WORKDIR/splitpipe" "$WORKDIR/rhapsody"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

ln -s "$PROJECT_ROOT/lib" "$WORKDIR/lib"
cat >"$WORKDIR/main.nf" <<'NF'
def here    = projectDir.toString()
def sp      = "${here}/splitpipe".toString()
def bd      = "${here}/rhapsody".toString()
def missing = "${here}/not_there".toString()

def base = [
    protocol: '10xv3', run_method: 'standard', perform_cellranger: false,
    splitpipe_installation: null, splitpipe_conda_env: null, splitpipe_demultiplex_script: null,
    splitpipe_kit: null, rhapsody_installation: null,
]
def ext = 'external_pipeline_only'

// [case, label, params, error substrings, warning substrings]; empty lists mean none
def cases = [
    ['external_only', '10x without Cell Ranger',   [protocol: '10xv3', run_method: ext],   ["Set 'perform_cellranger = true'"], []],
    ['external_only', 'OAK without Cell Ranger',   [protocol: 'oak_v1', run_method: ext],  ["Set 'perform_cellranger = true'"], []],
    ['external_only', 'Parse without split-pipe',  [protocol: 'parse_biosciences_WT', run_method: ext], ["Set 'splitpipe_installation'"], []],
    ['external_only', 'BD without its pipeline',   [protocol: 'bd_rhapsody_v1', run_method: ext], ["Set 'rhapsody_installation'"], []],
    ['external_only', 'MARS-seq',     [protocol: 'marsseq_v1', run_method: ext],           ["has no vendor pipeline"], []],
    ['external_only', 'sciRNAseq3',   [protocol: 'sciRNAseq3', run_method: ext],           ["has no vendor pipeline"], []],
    ['external_only', 'ScaleBio',     [protocol: 'scalebio_quantum_rna', run_method: ext], ["has no vendor pipeline"], []],
    ['external_only', 'seqspec',      [protocol: 'seqspec', run_method: ext],              ["has no vendor pipeline"], []],

    // A missing path is reported once, not also as "nothing to run"
    ['paths', 'Parse path missing (external)', [protocol: 'parse_biosciences_WT', run_method: ext, splitpipe_installation: missing, splitpipe_conda_env: 'env'],
        ["'splitpipe_installation' = '${missing}' does not exist"], []],
    ['paths', 'Parse path missing (standard)', [protocol: 'parse_biosciences_WT_mini', splitpipe_installation: missing, splitpipe_conda_env: 'env'],
        ["'splitpipe_installation' = '${missing}' does not exist"], []],
    ['paths', 'Parse without conda env', [protocol: 'parse_biosciences_WT', run_method: ext, splitpipe_installation: sp],
        ["'splitpipe_conda_env' is not"], []],
    ['paths', 'BD path missing (external)', [protocol: 'bd_rhapsody_enhancedbeads', run_method: ext, rhapsody_installation: missing],
        ["'rhapsody_installation' = '${missing}' does not exist"], []],
    ['paths', 'BD path missing (geneext_only)', [protocol: 'bd_rhapsody_v1', run_method: 'geneext_only', rhapsody_installation: missing],
        ["'rhapsody_installation' = '${missing}' does not exist"], []],

    ['mismatch', '10x with splitpipe_kit', [protocol: '10xv3', perform_cellranger: true, splitpipe_kit: 'WT'],
        [], ["'splitpipe_kit' is ignored: it belongs to split-pipe,\nwhich runs for parse_biosciences_WT_mini, parse_biosciences_WT, not for protocol '10xv3'"], 'docker'],
    ['mismatch', 'MARS-seq with perform_cellranger', [protocol: 'marsseq_v2', perform_cellranger: true],
        [], ["'perform_cellranger' is ignored"]],
    // Not this protocol's pipeline, so its path is not checked
    ['mismatch', 'Parse with a missing rhapsody_installation', [protocol: 'parse_biosciences_WT', rhapsody_installation: missing],
        [], ["'rhapsody_installation' is ignored"]],
    // The mismatch warns, and the run still stops for having nothing to run
    ['mismatch', 'MARS-seq external with perform_cellranger', [protocol: 'marsseq_v1', run_method: ext, perform_cellranger: true],
        ["has no vendor pipeline"], ["'perform_cellranger' is ignored"]],

    ['post_mapping', 'Cell Ranger and a missing split-pipe path', [protocol: 'parse_biosciences_WT', run_method: 'post_mapping', perform_cellranger: true, splitpipe_installation: missing],
        [], ["'perform_cellranger', 'splitpipe_installation' are ignored with 'run_method' = 'post_mapping'"]],
    ['post_mapping', 'nothing set', [protocol: 'sciRNAseq3', run_method: 'post_mapping'], [], []],

    // Each vendor pipeline under a profile it cannot run under; the sixth element is
    // workflow.containerEngine (null, the default, when no container engine is enabled)
    ['engine', 'Cell Ranger, no engine',  [protocol: '10xv3', perform_cellranger: true], [], ["Cell Ranger has no conda environment"]],
    ['engine', 'Cell Ranger, docker',     [protocol: '10xv3', perform_cellranger: true], [], [], 'docker'],
    ['engine', 'Cell Ranger off, no engine', [protocol: 'oak_v1'], [], []],
    ['engine', 'split-pipe, docker',      [protocol: 'parse_biosciences_WT', splitpipe_installation: sp, splitpipe_conda_env: 'env'],
        [], ["split-pipe has no container: under docker"], 'docker'],
    ['engine', 'split-pipe, no engine',   [protocol: 'parse_biosciences_WT', splitpipe_installation: sp, splitpipe_conda_env: 'env'], [], []],
    ['engine', 'BD, apptainer',           [protocol: 'bd_rhapsody_enhancedbeads', rhapsody_installation: bd],
        [], ["'cwltool --singularity', but under apptainer"], 'apptainer'],
    ['engine', 'BD, no engine',           [protocol: 'bd_rhapsody_v1', rhapsody_installation: bd], [], []],
    ['engine', 'BD off, singularity',     [protocol: 'bd_rhapsody_v1'], [], [], 'singularity'],
    // Vendor pipelines do not run with post_mapping, so only that warning is given
    ['engine', 'split-pipe post_mapping, docker', [protocol: 'parse_biosciences_WT', run_method: 'post_mapping', splitpipe_installation: sp, splitpipe_conda_env: 'env'],
        [], ["'splitpipe_installation', 'splitpipe_conda_env' are ignored with 'run_method' = 'post_mapping'"], 'docker'],

    ['clean', 'defaults',                 [:], [], []],
    ['clean', '10x external, Cell Ranger', [protocol: '10xv4', run_method: ext, perform_cellranger: true], [], [], 'singularity'],
    ['clean', 'Parse external, complete', [protocol: 'parse_biosciences_WT', run_method: ext, splitpipe_installation: sp, splitpipe_conda_env: 'env', splitpipe_kit: 'WT'], [], []],
    ['clean', 'BD external, complete',    [protocol: 'bd_rhapsody_v1', run_method: ext, rhapsody_installation: bd], [], []],
    ['clean', 'BD standard, nothing set', [protocol: 'bd_rhapsody_v1'], [], []],
    ['clean', 'MARS-seq standard',        [protocol: 'marsseq_v1'], [], []],
    ['clean', 'blank strings are unset',  [protocol: 'parse_biosciences_WT', splitpipe_installation: '  ', rhapsody_installation: ''], [], []],
    // Left to the error in preprocessing_workflow
    ['clean', 'unknown protocol',         [protocol: 'dropseq', run_method: ext], [], []],
]

cases.each { c ->
    def (name, label, overrides, want_errors, want_warnings) = c
    def engine = c.size() > 5 ? c[5] : null
    def got = VendorPipelines.check(base + overrides, engine)
    def problems = []
    [['error', got.errors, want_errors], ['warning', got.warnings, want_warnings]].each { kind, messages, wanted ->
        if (messages.size() != wanted.size()) {
            problems << "${messages.size()} ${kind}(s), want ${wanted.size()}: ${messages}"
        }
        wanted.each { text ->
            if (!messages.any { message -> message.contains(text.toString()) }) {
                problems << "no ${kind} with \"${text}\""
            }
        }
    }
    def detail = problems ? problems.join('; ') : (got.errors + got.warnings).collect { m -> m.readLines()[0] }.join(' | ')
    println "${problems ? 'BAD' : 'OK'}\t${name}\t${label}\t${detail}"
}
NF

LOG="$BCA_TEST_LOGDIR/vendor_pipelines.log"
if (cd "$WORKDIR" && run_logged "$LOG" nextflow -q run main.nf); then
    for case_name in external_only paths mismatch post_mapping engine clean; do
        expected=$(grep -c "^\s*\['$case_name'," "$WORKDIR/main.nf" || true)
        ok=$(grep -c -P "^OK\t$case_name\t" "$LOG" || true)
        bad=$(grep -P "^BAD\t$case_name\t" "$LOG" | head -2 | tr '\n' ' ' || true)
        if [[ -n "$bad" ]]; then
            record FAIL "$case_name" "$bad"
        elif [[ "$ok" -eq "$expected" ]]; then
            record PASS "$case_name" "$ok of $expected as expected"
        else
            record FAIL "$case_name" "expected $expected results, got $ok; see vendor_pipelines.log"
        fi
    done
else
    record FAIL "run" "nextflow run failed, see vendor_pipelines.log"
fi

finish_check
