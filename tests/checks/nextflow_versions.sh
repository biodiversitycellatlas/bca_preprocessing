#!/usr/bin/env bash
# description: Dry-run the pipeline against each installed Nextflow version
#
# Nextflow keeps every version it has been asked for under $NXF_HOME/framework,
# and NXF_VER decides which of them the launcher uses. This check walks those
# versions and puts the pipeline through increasingly demanding dry runs, so a
# version that no longer reads the config, load the plugins or resolve the
# workflow shows up here rather than at the start of a real run.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

NXF_FRAMEWORK_DIR="${NXF_HOME:-$HOME/.nextflow}/framework"

VERSIONS=""
DEPTH="parse"
PROFILE="conda"
CONFIG=""
INSTALL=0
IGNORE_MANIFEST=0
KEEP=0
TIMEOUT=600
LIST_ONLY=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh nextflow_versions [-- OPTIONS]
       tests/checks/nextflow_versions.sh [OPTIONS]

Check that the pipeline still works with each Nextflow version installed on this
machine. Nothing is executed, so no sequencing data is needed.

Options:
  --versions LIST     Comma-separated versions to test, e.g. 24.10.5,25.04.6
                      (default: every version under \$NXF_HOME/framework, plus
                      whatever \`nextflow -version\` reports)
  --install           Let Nextflow download versions that are not installed yet.
                      Without this, a version asked for by --versions but not
                      present is a failure.
  --depth LEVEL       How far to take each version. Each level includes the
                      previous ones:
                        launch   the version starts up at all
                        config   \`nextflow config\` resolves every config file
                        parse    plugins load and the entry script compiles
                                 (default)
                        preview  full dry run; resolves the channels of the
                                 whole workflow. Needs a config with real data
                                 paths, see --config.
  --profile NAME      Profile string passed to -profile (default: $PROFILE).
                      Use the same form as a real run, e.g. --profile crg,conda.
  --config FILE       Config used by the preview depth (default: the first
                      tests/*.config found).
  --ignore-manifest   Also test versions below the minimum declared by
                      manifest.nextflowVersion, which are skipped otherwise.
  --timeout SECONDS   Abort a single step after this long (default: $TIMEOUT,
                      0 disables). Needs the 'timeout' command.
  --keep              Keep the per-version working directories.
  --list-versions     Print the versions that would be tested and exit.
  -h, --help          Show this help.

Installing another version to test against:

  NXF_VER=24.10.5 nextflow -version      # downloads it and caches it
  tests/run_tests.sh nextflow_versions   # it is picked up from then on

Examples:
  # Everything installed, up to the entry script compiling
  tests/run_tests.sh nextflow_versions

  # Quick sweep: does each version still read the config?
  tests/run_tests.sh nextflow_versions -- --depth config

  # Full dry run of two specific versions, downloading them if needed
  tests/run_tests.sh nextflow_versions -- \\
      --versions 24.10.5,25.04.6 --install \\
      --depth preview --config tests/test_parsebio.config
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --versions)        VERSIONS="$2"; shift 2 ;;
        --install)         INSTALL=1; shift ;;
        --depth)           DEPTH="$2"; shift 2 ;;
        --profile)         PROFILE="$2"; shift 2 ;;
        --config)          CONFIG="$2"; shift 2 ;;
        --ignore-manifest) IGNORE_MANIFEST=1; shift ;;
        --timeout)         TIMEOUT="$2"; shift 2 ;;
        --keep)            KEEP=1; shift ;;
        --list-versions)   LIST_ONLY=1; shift ;;
        -h|--help)         usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage >&2; exit 2 ;;
    esac
done

# depth_rank LEVEL -- ordering for the cumulative --depth levels; -1 if unknown.
depth_rank() {
    case "$1" in
        launch)  printf '0' ;;
        config)  printf '1' ;;
        parse)   printf '2' ;;
        preview) printf '3' ;;
        *)       printf -- '-1' ;;
    esac
}

max_rank="$(depth_rank "$DEPTH")"
if [[ "$max_rank" -lt 0 ]]; then
    log_error "Invalid --depth: $DEPTH (expected launch, config, parse or preview)"
    exit 2
fi

# --------------------------------------------------------------------------
# Version discovery
# --------------------------------------------------------------------------

# installed_versions -- the versions cached under $NXF_HOME/framework.
installed_versions() {
    [[ -d "$NXF_FRAMEWORK_DIR" ]] || return 0
    find "$NXF_FRAMEWORK_DIR" -mindepth 1 -maxdepth 1 -type d | while IFS= read -r d; do
        basename "$d"
    done
}

# launcher_version -- the version the launcher uses when NXF_VER is unset.
# Never fails: it is called before we know whether nextflow is installed.
launcher_version() {
    have_cmd nextflow || return 0
    nextflow -version 2>/dev/null | sed -n 's/.*version \([0-9][^ ]*\).*/\1/p' | head -1 || true
}

# manifest_min -- the lowest version manifest.nextflowVersion allows, if it is
# expressed as a plain >= constraint.
manifest_min() {
    sed -n "s/.*nextflowVersion[^=]*=[^>]*>=[[:space:]]*\([0-9][0-9.]*\).*/\1/p" \
        "$PROJECT_ROOT/nextflow.config" | head -1
}

# version_lt A B -- true when A sorts before B.
version_lt() {
    [[ "$1" != "$2" && "$(printf '%s\n%s\n' "$1" "$2" | sort -V | head -1)" == "$1" ]]
}

LAUNCHER_VERSION="$(launcher_version)"

# version_available VERSION -- true when the launcher can use it without going
# to the network. A version installed as a self-contained distribution has no
# framework/ directory of its own, hence the second test.
version_available() {
    [[ -d "$NXF_FRAMEWORK_DIR/$1" || "$1" == "$LAUNCHER_VERSION" ]]
}

if [[ -n "$VERSIONS" ]]; then
    versions="$(split_csv "$VERSIONS")"
else
    versions="$({ installed_versions; printf '%s\n' "$LAUNCHER_VERSION"; } | sed '/^$/d' | sort -Vu)"
fi

if [[ "$LIST_ONLY" -eq 1 ]]; then
    if [[ -z "$versions" ]]; then
        log_warn "No Nextflow versions found under $NXF_FRAMEWORK_DIR"
    else
        printf '%s\n' "$versions"
    fi
    exit 0
fi

if ! have_cmd nextflow; then
    if [[ -n "$VERSIONS" ]]; then
        record FAIL "nextflow" "launcher not found on PATH"
    else
        record SKIP "nextflow" "launcher not installed"
    fi
    finish_check; exit $?
fi

if [[ -z "$versions" ]]; then
    record SKIP "nextflow" "no versions found under $NXF_FRAMEWORK_DIR"
    finish_check; exit $?
fi

min_version="$(manifest_min)"

# The preview depth is the only one that needs real data paths; it comes from
# the contributor's own config, since tests/*.config is not tracked in git.
preview_config="$CONFIG"
if [[ -z "$preview_config" && "$max_rank" -ge 3 ]]; then
    preview_config="$(find "$TESTS_DIR" -maxdepth 1 -name '*.config' | sort | head -1)"
fi
if [[ -n "$preview_config" && ! -f "$preview_config" ]]; then
    log_error "No such config file: $preview_config"
    exit 2
fi

n_versions="$(printf '%s\n' "$versions" | wc -l | tr -d ' ')"
log_step "Testing $n_versions Nextflow version(s) up to depth '$DEPTH' with -profile $PROFILE"
[[ -n "$min_version" ]] && log_dim "  manifest requires >=$min_version"

# --------------------------------------------------------------------------
# Per-version checks
# --------------------------------------------------------------------------

RUNDIR=""

# nxf VERSION LOGFILE ARGS... -- run the launcher pinned to VERSION, from the
# version's own working directory so .nextflow.log and work/ stay out of the
# repository.
nxf() {
    local version="$1" logfile="$2"; shift 2
    local -a cmd=(env "NXF_VER=$version" "NXF_ANSI_LOG=false" nextflow "$@")
    if [[ "$TIMEOUT" -gt 0 ]] && have_cmd timeout; then
        cmd=(timeout "$TIMEOUT" "${cmd[@]}")
    fi
    ( cd "$RUNDIR" && run_logged "$logfile" "${cmd[@]}" )
}

# step VERSION NAME ARGS... -- run one dry-run step and record the outcome.
# Returns non-zero when the step failed, so the caller can stop early.
step() {
    local version="$1" name="$2"; shift 2
    local logfile
    logfile="$BCA_TEST_LOGDIR/nextflow_versions__$(slugify "$version-$name").log"
    if nxf "$version" "$logfile" "$@"; then
        record PASS "$version: $name"
        return 0
    fi
    # Nextflow writes the interesting part of a failure to .nextflow.log, which
    # lives in the working directory rather than on stdout.
    [[ -f "$RUNDIR/.nextflow.log" ]] && cat "$RUNDIR/.nextflow.log" >>"$logfile"
    record FAIL "$version: $name" "see ${logfile#"$PROJECT_ROOT"/}"
    tail_log "$logfile"
    return 1
}

while IFS= read -r version; do
    [[ -n "$version" ]] || continue
    log_header "Nextflow $version"

    # Only reachable for a version named by --versions: discovered ones are, by
    # definition, already here.
    if [[ "$INSTALL" -eq 0 ]] && ! version_available "$version"; then
        record FAIL "$version" "not installed; add --install to download it"
        continue
    fi

    if [[ -n "$min_version" && "$IGNORE_MANIFEST" -eq 0 ]] && version_lt "$version" "$min_version"; then
        record SKIP "$version" "below manifest.nextflowVersion (>=$min_version)"
        continue
    fi

    RUNDIR="$BCA_TEST_LOGDIR/run_$(slugify "$version")"
    mkdir -p "$RUNDIR"

    step "$version" "launch" -version || continue

    [[ "$max_rank" -ge 1 ]] || continue
    step "$version" "config" config -profile "$PROFILE" "$PROJECT_ROOT" || continue

    # --help exits before any process runs, but only after the config is
    # applied, the plugins are resolved and the entry script is compiled.
    [[ "$max_rank" -ge 2 ]] || continue
    step "$version" "parse" run "$PROJECT_ROOT/main.nf" -profile "$PROFILE" --help || continue

    [[ "$max_rank" -ge 3 ]] || continue
    if [[ -z "$preview_config" ]]; then
        record SKIP "$version: preview" "no tests/*.config found; pass --config FILE"
        continue
    fi
    step "$version" "preview" run "$PROJECT_ROOT/main.nf" \
        -preview -profile "$PROFILE" -c "$preview_config" || continue
done <<<"$versions"

if [[ "$KEEP" -eq 0 ]]; then
    find "$BCA_TEST_LOGDIR" -maxdepth 1 -type d -name 'run_*' -exec rm -rf {} + 2>/dev/null || true
fi

finish_check
