#!/usr/bin/env bash
# description: Verify every module container image resolves for each container profile
#
# Extracts the container references declared by the modules and checks that each
# one can actually be resolved by the engines behind the container profiles in
# nextflow.config. This catches images that were deleted, retagged, or that a
# given engine simply cannot consume. It also checks that every module with a
# conda environment declares a container, and that each image provides what the
# module's environment.yml lists, so the container profiles run what conda runs.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

ALL_PROFILES="docker podman singularity apptainer wave"

PROFILES=""
DEPTH="auto"
ONLY=""
FORCE=0
CACHEDIR="$TESTS_DIR/.cache/containers"

usage() {
    cat <<EOF
Usage: tests/run_tests.sh containers [-- OPTIONS]
       tests/checks/containers.sh [OPTIONS]

Verify the container images declared by the modules against each container
profile in nextflow.config.

Options:
  --profiles LIST   Comma-separated subset of: $ALL_PROFILES
                    (default: every engine found on PATH)
  --depth LEVEL     How hard to check each image:
                      manifest  query the registry only; no layers downloaded
                      pull      download the full image
                      env       pull, then check that the image provides every
                                package in the module's environment.yml
                                (singularity/apptainer; needs python3)
                      auto      manifest for docker/podman, env for
                                singularity/apptainer (default)
  --only PATTERN    Only check images matching PATTERN (grep -E).
  --force           Re-pull images already present in the cache.
  --cachedir DIR    Where pulled images are stored (default: tests/.cache/containers).
  -h, --help        Show this help.

Notes:
  * The 'coverage' case needs no engine: it fails for every module that has an
    environment.yml but no container, since such a module cannot run under the
    container profiles.
  * Most modules declare oras:// images, which only singularity and apptainer
    can consume. They are reported as skipped under docker/podman rather than
    failed, with a count in the summary.
  * The wave profile builds images on demand through Seqera and needs
    TOWER_ACCESS_TOKEN plus the 'wave' CLI; without those it is skipped.

Examples:
  # Fast registry-only sweep across every available engine
  tests/run_tests.sh containers -- --depth manifest

  # Actually download the singularity images (also warms the cache)
  tests/run_tests.sh containers -- --profiles singularity --depth pull
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --profiles) PROFILES="$2"; shift 2 ;;
        --depth)    DEPTH="$2"; shift 2 ;;
        --only)     ONLY="$2"; shift 2 ;;
        --force)    FORCE=1; shift ;;
        --cachedir) CACHEDIR="$2"; shift 2 ;;
        -h|--help)  usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage >&2; exit 2 ;;
    esac
done

case "$DEPTH" in
    auto|manifest|pull|env) ;;
    *) log_error "Invalid --depth: $DEPTH (expected auto, manifest, pull or env)"; exit 2 ;;
esac

# --------------------------------------------------------------------------
# Coverage: a module with a conda environment needs a container as well
# --------------------------------------------------------------------------

# Without one, the task runs on the host under the container profiles and finds
# none of its tools. Needs no engine, so it runs on every machine.
log_header "coverage"
missing=""
while IFS= read -r yml; do
    nf="$(dirname "$yml")/main.nf"
    [[ -f "$nf" ]] || continue
    if ! grep -qE '^[[:space:]]*container[[:space:]]*["'"'"'{]' "$nf"; then
        module="${nf#"$PROJECT_ROOT"/modules/}"
        missing="$missing ${module%/main.nf}"
    fi
done < <(find "$PROJECT_ROOT/modules" -name environment.yml | sort)
if [[ -z "$missing" ]]; then
    record PASS "coverage" "every module with an environment.yml declares a container"
else
    record FAIL "coverage" "environment.yml without a container:$missing"
fi

# --------------------------------------------------------------------------
# Select engines
# --------------------------------------------------------------------------

if [[ -n "$PROFILES" ]]; then
    requested="$(split_csv "$PROFILES")"
    for p in $requested; do
        if ! printf '%s\n' $ALL_PROFILES | grep -qx "$p"; then
            log_error "Not a supported container profile: $p (expected one of: $ALL_PROFILES)"
            exit 2
        fi
    done
else
    requested="$ALL_PROFILES"
fi

engines=""
for p in $requested; do
    case "$p" in
        wave)
            if have_cmd wave && [[ -n "${TOWER_ACCESS_TOKEN:-}" ]]; then
                engines="$engines wave"
            else
                record SKIP "wave" "needs the 'wave' CLI and TOWER_ACCESS_TOKEN"
            fi
            ;;
        *)
            if have_cmd "$p"; then
                engines="$engines $p"
            elif [[ -n "$PROFILES" ]]; then
                record FAIL "$p" "engine not found on PATH"
            else
                record SKIP "$p" "engine not installed"
            fi
            ;;
    esac
done

if [[ -z "${engines// /}" ]]; then
    log_warn "No container engines available; nothing to do."
    finish_check; exit $?
fi

# --------------------------------------------------------------------------
# Collect container references
# --------------------------------------------------------------------------

refs_file="$(mktemp -t bca_refs.XXXXXX)"
trap 'rm -f "$refs_file"' EXIT

n_dynamic=0
while IFS= read -r nf; do
    grep -qE '^[[:space:]]*container[[:space:]]*["'"'"'{]' "$nf" || continue
    module="${nf#"$PROJECT_ROOT"/modules/}"; module="${module%/main.nf}"
    literals="$(awk -f "$BCA_LIB_DIR/extract_containers.awk" "$nf")"
    if [[ -z "$literals" ]]; then
        # e.g. `container { demuxafy_sif }` -- resolved at runtime from params,
        # so there is no static reference to validate.
        record SKIP "$module" "container resolved at runtime from params"
        n_dynamic=$((n_dynamic + 1))
        continue
    fi
    while IFS= read -r ref; do
        # Redirect per line, not on the enclosing loop: the loop's stdout has to
        # stay free for record() to report on.
        if [[ -n "$ref" ]]; then
            printf '%s\t%s\n' "$ref" "$module" >>"$refs_file"
        fi
    done <<<"$literals"
done < <(find "$PROJECT_ROOT/modules" -name main.nf | sort)

# Several modules share an image (four samtools modules use one tag): each unique
# reference is pulled once, and checked against the environment of every module
# that declares it.
sort -u -o "$refs_file" "$refs_file"

if [[ -n "$ONLY" ]]; then
    grep -E "$ONLY" "$refs_file" >"$refs_file.tmp" || true
    mv "$refs_file.tmp" "$refs_file"
fi

n_refs="$(cut -f1 "$refs_file" | sort -u | wc -l | tr -d ' ')"
if [[ "$n_refs" -eq 0 ]]; then
    log_error "No container references matched${ONLY:+ pattern '$ONLY'}."
    exit 1
fi

log_step "Checking $n_refs unique image(s) with:${engines}"
[[ "$n_dynamic" -gt 0 ]] && log_dim "  $n_dynamic module(s) use a runtime-resolved container and were skipped"

mkdir -p "$CACHEDIR"

# --------------------------------------------------------------------------
# Per-engine checks
# --------------------------------------------------------------------------

# effective_depth ENGINE -- resolve "auto" into a concrete depth.
effective_depth() {
    if [[ "$DEPTH" != auto ]]; then printf '%s' "$DEPTH"; return; fi
    case "$1" in
        docker|podman) printf 'manifest' ;;
        *)             printf 'env' ;;
    esac
}

# check_env ENGINE SIF REF MODULES -- check a pulled image against the
# environment.yml of each module that declares it (see lib/check_image_env.py).
check_env() {
    local engine="$1" sif="$2" ref="$3" modules="$4"
    local module yml listing logfile
    if ! have_cmd python3; then
        record SKIP "$engine env: $ref" "needs python3"
        return
    fi
    listing="$(mktemp -t bca_listing.XXXXXX)"
    # A glob that matches nothing makes ls fail, so the listing is judged by its content
    "$engine" exec "$sif" sh -c 'ls -1 /opt/conda/conda-meta/*.json /opt/conda/lib/python3*/site-packages/ 2>/dev/null' \
        >"$listing" 2>/dev/null || true
    for module in $modules; do
        yml="$PROJECT_ROOT/modules/$module/environment.yml"
        [[ -f "$yml" ]] || continue
        logfile="$BCA_TEST_LOGDIR/containers__env_$(slugify "$engine-$module").log"
        if run_logged "$logfile" python3 "$BCA_LIB_DIR/check_image_env.py" "$yml" "$listing"; then
            record PASS "$engine env: $module" "$(tail -n1 "$logfile")"
        else
            record FAIL "$engine env: $module" "$(grep -cE '^(missing|mismatch)' "$logfile") problem(s) in ${ref##*/}"
            grep -E '^(missing|mismatch)|^No conda' "$logfile" | head -n 12 | sed 's/^/        | /'
        fi
    done
    rm -f "$listing"
}

# check_ref ENGINE REF MODULES -- MODULES lists every module declaring REF
check_ref() {
    local engine="$1" ref="$2" module="$3"
    local depth logfile slug sif
    depth="$(effective_depth "$engine")"
    slug="$(slugify "$engine-$ref")"
    logfile="$BCA_TEST_LOGDIR/containers__$slug.log"

    case "$engine" in
        docker|podman)
            case "$ref" in
                oras://*)
                    record SKIP "$engine: $ref" "oras:// artifact, singularity/apptainer only"
                    return ;;
                https://*|http://*)
                    # A raw singularity image URL; nothing for docker to do.
                    record SKIP "$engine: $ref" "direct image URL, singularity/apptainer only"
                    return ;;
            esac
            local image="${ref#docker://}"
            if [[ "$depth" == manifest ]]; then
                if run_logged "$logfile" "$engine" manifest inspect "$image"; then
                    record PASS "$engine: $image"
                else
                    record FAIL "$engine: $image" "manifest not resolvable ($module)"
                    tail_log "$logfile"
                fi
            else
                if run_logged "$logfile" "$engine" pull "$image"; then
                    record PASS "$engine: $image"
                else
                    record FAIL "$engine: $image" "pull failed ($module)"
                    tail_log "$logfile"
                fi
            fi
            ;;

        singularity|apptainer)
            case "$ref" in
                https://*|http://*)
                    if run_logged "$logfile" curl -fsSL -I "$ref"; then
                        record PASS "$engine: $ref"
                    else
                        record FAIL "$engine: $ref" "URL not reachable ($module)"
                        tail_log "$logfile"
                    fi
                    return ;;
            esac

            # singularity has no way to inspect a remote image without fetching
            # it, so a manifest-depth run borrows docker for registry images and
            # skips oras:// ones.
            if [[ "$depth" == manifest ]]; then
                case "$ref" in
                    oras://*)
                        record SKIP "$engine: $ref" "oras:// cannot be checked without pulling; use --depth pull"
                        return ;;
                esac
                if have_cmd docker; then
                    if run_logged "$logfile" docker manifest inspect "${ref#docker://}"; then
                        record PASS "$engine: $ref" "via docker manifest"
                    else
                        record FAIL "$engine: $ref" "manifest not resolvable ($module)"
                        tail_log "$logfile"
                    fi
                else
                    record SKIP "$engine: $ref" "manifest depth needs docker; use --depth pull"
                fi
                return
            fi

            sif="$CACHEDIR/$(slugify "$ref").sif"
            if [[ -s "$sif" && "$FORCE" -eq 0 ]]; then
                record PASS "$engine: $ref" "cached"
            # A plain registry reference needs a transport; Nextflow adds the same one
            elif run_logged "$logfile" "$engine" pull --force "$sif" "$([[ "$ref" == *://* ]] || printf 'docker://')$ref"; then
                record PASS "$engine: $ref" "$(du -h "$sif" 2>/dev/null | cut -f1)"
            else
                record FAIL "$engine: $ref" "pull failed ($module)"
                rm -f "$sif"
                tail_log "$logfile"
                return
            fi
            if [[ "$depth" == env ]]; then
                check_env "$engine" "$sif" "$ref" "$module"
            fi
            ;;

        wave)
            # Ask Wave to resolve/build the image; this is what the wave profile
            # does at run time, minus the pipeline execution.
            if run_logged "$logfile" wave -i "${ref#oras://}"; then
                record PASS "wave: $ref"
            else
                record FAIL "wave: $ref" "wave could not resolve the image ($module)"
                tail_log "$logfile"
            fi
            ;;
    esac
}

for engine in $engines; do
    log_header "$engine (depth: $(effective_depth "$engine"))"
    while IFS= read -r ref; do
        [[ -n "$ref" ]] || continue
        modules="$(awk -F'\t' -v r="$ref" '$1 == r { printf "%s ", $2 }' "$refs_file")"
        check_ref "$engine" "$ref" "${modules% }"
    done < <(cut -f1 "$refs_file" | sort -u)
done

finish_check
