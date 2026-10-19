#!/usr/bin/env bash
# Entry point for the BCA pre-processing test suite.
#
# Every executable in tests/checks/ is a check and is discovered automatically,
# so adding a new one is a matter of dropping in a script -- see tests/README.md.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/lib" && pwd)/common.sh"

CHECKS_DIR="$TESTS_DIR/checks"

usage() {
    cat <<EOF
Usage: tests/run_tests.sh [CHECK...] [-- OPTIONS]

Run the pipeline test suite. With no arguments, every check runs with its
default settings.

Arguments:
  CHECK...      Names of checks to run (default: all). See --list.
  -- OPTIONS    Everything after -- is passed through to each selected check.

Options:
  -l, --list    List the available checks and exit.
  -h, --help    Show this help.

Examples:
  tests/run_tests.sh                                  # everything
  tests/run_tests.sh --list
  tests/run_tests.sh conda_envs                       # one check
  tests/run_tests.sh conda_envs -- --solve-only       # with check options
  tests/run_tests.sh containers -- --depth manifest

Logs for each run are written to tests/.logs/<timestamp>/.
EOF
}

# discover -- print the name of every check, one per line, in run order.
discover() {
    find "$CHECKS_DIR" -maxdepth 1 -name '*.sh' | sort | while IFS= read -r f; do
        basename "$f" .sh
    done
}

# describe CHECK -- the check's `# description:` line.
describe() {
    sed -n 's/^# description: *//p' "$CHECKS_DIR/$1.sh" | head -1
}

list_checks() {
    printf '%sAvailable checks%s\n\n' "$C_BOLD" "$C_RESET"
    while IFS= read -r name; do
        printf '  %-18s %s\n' "$name" "$(describe "$name")"
    done < <(discover)
    printf '\nRun a check with:  tests/run_tests.sh <check>\n'
    printf 'See its options:   tests/checks/<check>.sh --help\n'
}

# --------------------------------------------------------------------------
# Parse arguments
# --------------------------------------------------------------------------

selected=()
passthrough=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -l|--list) list_checks; exit 0 ;;
        -h|--help) usage; exit 0 ;;
        --) shift; passthrough=("$@"); break ;;
        -*) log_error "Unknown option: $1"; usage >&2; exit 2 ;;
        *)  selected+=("$1"); shift ;;
    esac
done

available="$(discover)"
if [[ ${#selected[@]} -eq 0 ]]; then
    while IFS= read -r name; do selected+=("$name"); done <<<"$available"
else
    for name in "${selected[@]}"; do
        if ! printf '%s\n' "$available" | grep -qx "$name"; then
            log_error "No such check: $name"
            printf '\n' >&2; list_checks >&2
            exit 2
        fi
    done
fi

# --------------------------------------------------------------------------
# Run the selected checks
# --------------------------------------------------------------------------

: >"$BCA_TEST_RESULTS"

printf '%sBCA pre-processing test suite%s\n' "$C_BOLD" "$C_RESET"
printf 'Checks: %s\n' "${selected[*]}"
printf 'Logs:   %s\n' "${BCA_TEST_LOGDIR#"$PROJECT_ROOT"/}"

failed_checks=()
start_all="$SECONDS"

for name in "${selected[@]}"; do
    printf '\n%s%s %s %s%s\n' "$C_BOLD$C_BLUE" "========" "$name" "========" "$C_RESET"
    printf '%s%s%s\n' "$C_DIM" "$(describe "$name")" "$C_RESET"

    # Each check runs in its own process so one failing check cannot abort the
    # run; results reach us through the shared $BCA_TEST_RESULTS file.
    if BCA_CHECK_NAME="$name" bash "$CHECKS_DIR/$name.sh" "${passthrough[@]+"${passthrough[@]}"}"; then
        :
    else
        failed_checks+=("$name")
    fi
done

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------

printf '\n%s%s%s\n' "$C_BOLD" "=========== summary ===========" "$C_RESET"
printf '\n%-18s %7s %7s %7s\n' "CHECK" "PASS" "FAIL" "SKIP"
for name in "${selected[@]}"; do
    p="$(count_results PASS "$name")"
    f="$(count_results FAIL "$name")"
    s="$(count_results SKIP "$name")"
    if [[ "$f" -gt 0 ]]; then colour="$C_RED"; else colour="$C_GREEN"; fi
    printf '%-18s %s%7s%s %s%7s%s %7s\n' "$name" "$C_GREEN" "$p" "$C_RESET" "$colour" "$f" "$C_RESET" "$s"
done

total_fail="$(awk -F'\t' '$2=="FAIL"' "$BCA_TEST_RESULTS" | wc -l | tr -d ' ')"
printf '\nTook %s\n' "$(human_time $((SECONDS - start_all)))"

if [[ "$total_fail" -gt 0 ]]; then
    printf '\n%sFailures:%s\n' "$C_BOLD$C_RED" "$C_RESET"
    awk -F'\t' '$2=="FAIL" { printf "  %s: %s  %s\n", $1, $3, $4 }' "$BCA_TEST_RESULTS"
    printf '\nFull logs in %s\n' "${BCA_TEST_LOGDIR#"$PROJECT_ROOT"/}"
    exit 1
fi

if [[ ${#failed_checks[@]} -gt 0 ]]; then
    log_error "Checks reported an error: ${failed_checks[*]}"
    exit 1
fi

printf '\n%sAll checks passed.%s\n' "$C_GREEN$C_BOLD" "$C_RESET"
