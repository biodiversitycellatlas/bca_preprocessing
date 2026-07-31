#!/usr/bin/env bash
# Shared helpers for the BCA pre-processing test suite.
#
# Sourced by tests/run_tests.sh and by every check in tests/checks/. Checks are
# plain executables, so they can also be sourced standalone:
#     source "$(dirname "$0")/../lib/common.sh"

# Resolve locations from this file, so checks work regardless of the caller's cwd.
BCA_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTS_DIR="$(cd "$BCA_LIB_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$TESTS_DIR/.." && pwd)"
export BCA_LIB_DIR TESTS_DIR PROJECT_ROOT

# Every check appends its results here as TSV, so run_tests.sh can summarise
# across checks. Set by the runner; created on demand for standalone runs.
if [[ -z "${BCA_TEST_RESULTS:-}" ]]; then
    BCA_TEST_RESULTS="$(mktemp -t bca_test_results.XXXXXX)"
    export BCA_TEST_RESULTS
fi

# Per-run log directory. One log file per individual test case.
if [[ -z "${BCA_TEST_LOGDIR:-}" ]]; then
    BCA_TEST_LOGDIR="$TESTS_DIR/.logs/$(date +%Y%m%d_%H%M%S)"
    export BCA_TEST_LOGDIR
fi
mkdir -p "$BCA_TEST_LOGDIR"

# Name of the check currently running, used to group results in the summary.
BCA_CHECK_NAME="${BCA_CHECK_NAME:-$(basename "${0%.sh}")}"
export BCA_CHECK_NAME

# --------------------------------------------------------------------------
# Output formatting
# --------------------------------------------------------------------------

if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
    C_RESET=$'\033[0m'; C_RED=$'\033[31m'; C_GREEN=$'\033[32m'
    C_YELLOW=$'\033[33m'; C_BLUE=$'\033[34m'; C_BOLD=$'\033[1m'; C_DIM=$'\033[2m'
else
    C_RESET=''; C_RED=''; C_GREEN=''; C_YELLOW=''; C_BLUE=''; C_BOLD=''; C_DIM=''
fi

log_info()  { printf '%s\n' "$*"; }
log_step()  { printf '%s==>%s %s\n' "$C_BLUE$C_BOLD" "$C_RESET" "$*"; }
log_warn()  { printf '%sWARN%s %s\n' "$C_YELLOW" "$C_RESET" "$*" >&2; }
log_error() { printf '%sERROR%s %s\n' "$C_RED" "$C_RESET" "$*" >&2; }
log_dim()   { printf '%s%s%s\n' "$C_DIM" "$*" "$C_RESET"; }

log_header() {
    printf '\n%s%s%s\n' "$C_BOLD" "$*" "$C_RESET"
    printf '%s\n' "$(printf '%.0s-' $(seq 1 ${#1}))"
}

# --------------------------------------------------------------------------
# Result recording
# --------------------------------------------------------------------------

# record STATUS NAME [DETAIL]
#
# STATUS is PASS, FAIL or SKIP. Results are appended to $BCA_TEST_RESULTS so
# that they survive subshells -- checks may run cases in parallel, where
# in-memory counters in the parent shell would be lost.
record() {
    local status="$1" name="$2" detail="${3:-}"
    printf '%s\t%s\t%s\t%s\n' "$BCA_CHECK_NAME" "$status" "$name" "$detail" >>"$BCA_TEST_RESULTS"
    case "$status" in
        PASS) printf '  %sPASS%s  %s%s\n' "$C_GREEN" "$C_RESET" "$name" "${detail:+  $C_DIM$detail$C_RESET}" ;;
        FAIL) printf '  %sFAIL%s  %s%s\n' "$C_RED" "$C_RESET" "$name" "${detail:+  $detail}" ;;
        SKIP) printf '  %sSKIP%s  %s%s\n' "$C_YELLOW" "$C_RESET" "$name" "${detail:+  $C_DIM$detail$C_RESET}" ;;
    esac
}

# count_results STATUS [CHECK] -- number of recorded results with that status.
count_results() {
    local status="$1" check="${2:-$BCA_CHECK_NAME}"
    awk -F'\t' -v c="$check" -v s="$status" '$1==c && $2==s' "$BCA_TEST_RESULTS" | wc -l | tr -d ' '
}

# finish_check -- print the per-check tally and exit non-zero if anything failed.
finish_check() {
    local pass fail skip
    pass="$(count_results PASS)"; fail="$(count_results FAIL)"; skip="$(count_results SKIP)"
    printf '\n%s%s%s: %s%s passed%s, %s%s failed%s, %s%s skipped%s\n' \
        "$C_BOLD" "$BCA_CHECK_NAME" "$C_RESET" \
        "$C_GREEN" "$pass" "$C_RESET" \
        "$([[ $fail -gt 0 ]] && printf '%s' "$C_RED")" "$fail" "$C_RESET" \
        "$C_YELLOW" "$skip" "$C_RESET"
    [[ "$fail" -eq 0 ]]
}

# --------------------------------------------------------------------------
# Misc helpers
# --------------------------------------------------------------------------

have_cmd() { command -v "$1" >/dev/null 2>&1; }

# slugify STRING -- filesystem-safe identifier, used for log and cache names.
slugify() { printf '%s' "$1" | tr -c 'A-Za-z0-9._-' '_'; }

# human_time SECONDS -- e.g. "3m12s"
human_time() {
    local s="$1"
    if [[ "$s" -ge 60 ]]; then printf '%dm%02ds' $((s / 60)) $((s % 60)); else printf '%ds' "$s"; fi
}

# split_csv "a,b,c" -- print one item per line.
split_csv() { printf '%s' "$1" | tr ',' '\n' | sed '/^$/d'; }

# run_logged LOGFILE COMMAND... -- run a command, tee-ing all output to a log
# file. Returns the command's exit status.
run_logged() {
    local logfile="$1"; shift
    printf '$ %s\n\n' "$*" >"$logfile"
    "$@" >>"$logfile" 2>&1
}

# tail_log LOGFILE [N] -- last N lines of a log, indented, for failure context.
tail_log() {
    local logfile="$1" n="${2:-12}"
    [[ -f "$logfile" ]] || return 0
    printf '%s' "$C_DIM"
    tail -n "$n" "$logfile" | sed 's/^/        | /'
    printf '%s' "$C_RESET"
}
