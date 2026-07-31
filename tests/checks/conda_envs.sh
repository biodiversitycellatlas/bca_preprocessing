#!/usr/bin/env bash
# description: Build every module environment.yml with each conda-based profile
#
# Walks every environment.yml under modules/ and creates it with the solver
# behind each conda-based profile in nextflow.config (conda, mamba, micromamba).
# This catches unsolvable pins and packages that have been removed from a
# channel before they surface halfway through a pipeline run.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

ALL_PROFILES="conda mamba micromamba"

PROFILES=""
SOLVE_ONLY=0
KEEP=0
JOBS=1
ONLY=""
ENVDIR="$TESTS_DIR/.envs"

usage() {
    cat <<EOF
Usage: tests/run_tests.sh conda_envs [-- OPTIONS]
       tests/checks/conda_envs.sh [OPTIONS]

Create every modules/**/environment.yml with each conda-based profile.

Options:
  --profiles LIST   Comma-separated subset of: $ALL_PROFILES
                    (default: every solver found on PATH)
  --solve-only      Resolve dependencies without downloading or installing.
                    Much faster; catches unsatisfiable specs but not broken
                    packages. Needs a solver new enough to support --dry-run.
  --keep            Keep the created environments instead of removing each one
                    after it is verified. Needs a lot of disk space.
  --jobs N          Create N environments in parallel (default: $JOBS). Output
                    from parallel cases interleaves.
  --only PATTERN    Only test modules whose path matches PATTERN (grep -E).
  --envdir DIR      Where to create environments (default: tests/.envs).
  -h, --help        Show this help.

Examples:
  # Fast sanity check that every recipe still solves
  tests/run_tests.sh conda_envs -- --solve-only

  # Full build of just the samtools modules, under mamba only
  tests/run_tests.sh conda_envs -- --profiles mamba --only samtools
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --profiles)   PROFILES="$2"; shift 2 ;;
        --solve-only) SOLVE_ONLY=1; shift ;;
        --keep)       KEEP=1; shift ;;
        --jobs)       JOBS="$2"; shift 2 ;;
        --only)       ONLY="$2"; shift 2 ;;
        --envdir)     ENVDIR="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage >&2; exit 2 ;;
    esac
done

# --------------------------------------------------------------------------
# Select solvers
# --------------------------------------------------------------------------

if [[ -n "$PROFILES" ]]; then
    requested="$(split_csv "$PROFILES")"
    for p in $requested; do
        if ! printf '%s\n' $ALL_PROFILES | grep -qx "$p"; then
            log_error "Not a conda-based profile: $p (expected one of: $ALL_PROFILES)"
            exit 2
        fi
    done
else
    requested="$ALL_PROFILES"
fi

solvers=""
for p in $requested; do
    if have_cmd "$p"; then
        solvers="$solvers $p"
    elif [[ -n "$PROFILES" ]]; then
        # Explicitly asked for -- absence is a failure, not a silent skip.
        record FAIL "$p" "solver not found on PATH"
    else
        record SKIP "$p" "solver not installed"
    fi
done

if [[ -z "${solvers// /}" ]]; then
    log_warn "No conda solvers available; nothing to do."
    finish_check; exit $?
fi

# --------------------------------------------------------------------------
# Discover environment files
# --------------------------------------------------------------------------

env_files="$(find "$PROJECT_ROOT/modules" -name environment.yml | sort)"
if [[ -n "$ONLY" ]]; then
    env_files="$(printf '%s\n' "$env_files" | grep -E "$ONLY" || true)"
fi

if [[ -z "$env_files" ]]; then
    log_error "No environment.yml files matched${ONLY:+ pattern '$ONLY'}."
    exit 1
fi

n_envs="$(printf '%s\n' "$env_files" | wc -l | tr -d ' ')"
log_step "Testing $n_envs environment(s) with:${solvers}"
[[ "$SOLVE_ONLY" -eq 1 ]] && log_dim "  solve-only mode: dependencies are resolved but not installed"
[[ "$KEEP" -eq 1 ]] && log_dim "  keeping environments under $ENVDIR"

mkdir -p "$ENVDIR"

# --------------------------------------------------------------------------
# Build one environment
# --------------------------------------------------------------------------

build_env() {
    local solver="$1" yml="$2"
    local module slug prefix logfile start elapsed
    module="${yml#"$PROJECT_ROOT"/modules/}"; module="${module%/environment.yml}"
    slug="$(slugify "$solver-$module")"
    prefix="$ENVDIR/$slug"
    logfile="$BCA_TEST_LOGDIR/conda_envs__$slug.log"

    rm -rf "$prefix"

    local -a cmd
    case "$solver" in
        conda)      cmd=(conda env create --quiet --file "$yml" --prefix "$prefix") ;;
        mamba)      cmd=(mamba env create --quiet --file "$yml" --prefix "$prefix") ;;
        micromamba) cmd=(micromamba create --quiet --yes --file "$yml" --prefix "$prefix") ;;
    esac
    [[ "$SOLVE_ONLY" -eq 1 ]] && cmd+=(--dry-run)

    start="$SECONDS"
    if run_logged "$logfile" "${cmd[@]}"; then
        elapsed=$((SECONDS - start))
        record PASS "$solver: $module" "$(human_time "$elapsed")"
        [[ "$KEEP" -eq 1 || "$SOLVE_ONLY" -eq 1 ]] || rm -rf "$prefix"
    else
        elapsed=$((SECONDS - start))
        record FAIL "$solver: $module" "see ${logfile#"$PROJECT_ROOT"/}"
        tail_log "$logfile"
        rm -rf "$prefix"
    fi
}

# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------

for solver in $solvers; do
    log_header "$solver"
    running=0
    while IFS= read -r yml; do
        [[ -n "$yml" ]] || continue
        if [[ "$JOBS" -gt 1 ]]; then
            build_env "$solver" "$yml" &
            running=$((running + 1))
            if [[ "$running" -ge "$JOBS" ]]; then wait -n 2>/dev/null || wait; running=$((running - 1)); fi
        else
            build_env "$solver" "$yml"
        fi
    done <<<"$env_files"
    wait
done

[[ "$KEEP" -eq 1 ]] || rmdir "$ENVDIR" 2>/dev/null || true

finish_check
