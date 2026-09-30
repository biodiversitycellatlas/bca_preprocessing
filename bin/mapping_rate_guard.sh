#!/usr/bin/env bash
#
# Mapping-rate guard for STARsolo and alevin-fry (salmon alevin).
#
# A sample whose reads mostly fail to map (contamination, wrong species or reference,
# wrong protocol geometry) would otherwise spend hours in the mapper and then run every
# downstream step on almost nothing. This wraps the mapper command, watches its mapping
# rate and cancels the task with exit code 42 once the rate is known to be below
# the threshold. The errorStrategy in lib/MappingRateCheck.groovy turns that exit code
# into either a pipeline abort or a dropped sample.
#
#   STAR    Log.progress.out is refreshed once a minute while STAR maps. The rate is
#           checked once --check-reads reads have been processed, and again from
#           Log.final.out when STAR finishes if no verdict was reached before.
#   salmon  writes no progress file that can be parsed, so its rate is read from
#           aux_info/meta_info.json as soon as `salmon alevin --justAlign` finishes,
#           before alevin-fry's permit-list/collate/quant steps.
#
# When both mappers run on one sample they share a verdict file in --flag-dir, keyed
# on the sample's base id. The first verdict wins: a "fail" makes the other mapper kill
# itself (or never start), a "pass" lets it skip its own check.
#
# Only bash, awk, grep and sed are used: the STAR container has no python or jq.
#
# Usage:
#   mapping_rate_guard.sh run --mapper star|salmon --flag-dir DIR --sample BASE_ID \
#       --label TASK_ID --min PCT --report FILE [--progress FILE] [--check-reads N] \
#       [--poll SECS] -- COMMAND [ARGS...]
#   mapping_rate_guard.sh star-progress Log.progress.out    -> "reads pct"
#   mapping_rate_guard.sh star-final    Log.final.out       -> "reads pct"
#   mapping_rate_guard.sh salmon        meta_info.json      -> "reads pct"
#   mapping_rate_guard.sh verdict-get   DIR BASE_ID         -> pass | fail | none
#   mapping_rate_guard.sh verdict-set   DIR BASE_ID pass|fail MAPPER LABEL STAGE READS PCT MIN

set -uo pipefail

EXIT_LOW_MAPPING=42

# --------------------------------------------------------------------------
# Parsers. Each prints "reads pct", or nothing when the value is not there yet.
# --------------------------------------------------------------------------

# Data rows start with a timestamp ("Sep 30 10:01:23"), followed by speed, read
# number, read length and "Mapped unique" as the first percentage. Header rows and
# the closing "ALL DONE!" line have no HH:MM:SS field and are skipped.
star_progress() {
    [[ -s "$1" ]] || return 0
    awk '
        $3 ~ /^[0-9][0-9]:[0-9][0-9]:[0-9][0-9]$/ {
            for (i = 4; i <= NF; i++) if ($i ~ /%$/) {
                pct = $i; sub(/%$/, "", pct); reads = $(i - 2); break
            }
        }
        END { if (reads != "") print reads, pct }
    ' "$1"
}

star_final() {
    [[ -s "$1" ]] || return 0
    awk -F'|' '
        { key = $1; gsub(/^[ \t]+|[ \t]+$/, "", key); val = $2; gsub(/[ \t%]/, "", val) }
        key == "Number of input reads"   { reads = val }
        key == "Uniquely mapped reads %" { pct = val }
        END { if (reads != "" && pct != "") print reads, pct }
    ' "$1"
}

salmon_meta() {
    [[ -s "$1" ]] || return 0
    local reads pct
    reads=$(grep -o '"num_processed"[^,}]*' "$1" | head -n1 | sed 's/.*:[[:space:]]*//')
    pct=$(grep -o '"percent_mapped"[^,}]*' "$1" | head -n1 | sed 's/.*:[[:space:]]*//')
    [[ -n "$reads" && -n "$pct" ]] && echo "$reads $pct"
}

# --------------------------------------------------------------------------
# Shared verdict between the mappers of one sample
# --------------------------------------------------------------------------

verdict_get() {
    local file="$1/$2.verdict"
    if [[ -s "$file" ]]; then cut -f2 "$file"; else echo none; fi
}

# The first writer wins: the record is written to a private file and hard-linked into
# place, which fails atomically if another mapper got there first (also on NFS).
verdict_set() {
    local dir="$1" id="$2" verdict="$3" mapper="$4" label="$5" stage="$6" reads="$7" pct="$8" min="$9"
    mkdir -p "$dir"
    local tmp="$dir/.$id.$$.$RANDOM"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$id" "$verdict" "$mapper" "$label" "$stage" "$reads" "$pct" "$min" >"$tmp"
    ln "$tmp" "$dir/$id.verdict" 2>/dev/null
    local rc=$?
    rm -f "$tmp"
    return $rc
}

# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

mapper_name() { [[ "$1" == "star" ]] && echo "STARsolo" || echo "alevin-fry (salmon)"; }
rate_name()   { [[ "$1" == "star" ]] && echo "uniquely mapped" || echo "mapped"; }

# pct_below PCT MIN -- floating-point comparison
pct_below() { awk -v p="$1" -v m="$2" 'BEGIN { exit !(p + 0 < m + 0) }'; }

human_reads() { awk -v r="$1" 'BEGIN { if (r >= 1e6) printf "%.1fM", r / 1e6; else printf "%d", r }'; }

fail_banner() {
    local mapper="$1" label="$2" stage="$3" reads="$4" pct="$5" min="$6"
    {
        echo ""
        echo "==================== MAPPING RATE CHECK FAILED ===================="
        if [[ "$stage" == "other-mapper" ]]; then
            echo "Sample ${label}: cancelled because the other mapper already found this"
            echo "sample below the ${min}% mapping-rate threshold (see its task log and"
            echo "pipeline_info/mapping_rate_check.tsv)."
        else
            echo "Sample ${label} ($(mapper_name "$mapper"), ${stage} check after $(human_reads "$reads") reads):"
            echo "  ${pct}% of reads $(rate_name "$mapper"), below the ${min}% threshold."
        fi
        echo "Mapping and all downstream steps for this sample have been cancelled."
        echo ""
        echo "Likely causes: contamination, wrong reference/species, or a wrong"
        echo "protocol / barcode geometry for this library."
        echo ""
        echo "Recommendation: classify the unmapped reads with Kraken2 to find out what"
        echo "they are. Rerun this sample with the check disabled and Kraken enabled:"
        echo "   skip_mapping_rate_check = true"
        echo "   perform_kraken = true"
        echo "   star_generateBAM = true"
        echo "or run kraken2 directly on a subsample of the cDNA FASTQ."
        echo ""
        echo "To disable this check: skip_mapping_rate_check = true"
        echo "To lower the threshold: min_mapping_rate = <percent>"
        echo "==================================================================="
        echo ""
    } >&2
}

# --------------------------------------------------------------------------
# run: wrap the mapper
# --------------------------------------------------------------------------

run() {
    local mapper="" flag_dir="" sample="" label="" min="" report="" progress=""
    local check_reads=0 poll=60
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --mapper)      mapper="$2"; shift 2 ;;
            --flag-dir)    flag_dir="$2"; shift 2 ;;
            --sample)      sample="$2"; shift 2 ;;
            --label)       label="$2"; shift 2 ;;
            --min)         min="$2"; shift 2 ;;
            --report)      report="$2"; shift 2 ;;
            --progress)    progress="$2"; shift 2 ;;
            --check-reads) check_reads="$2"; shift 2 ;;
            --poll)        poll="$2"; shift 2 ;;
            --)            shift; break ;;
            *) echo "mapping_rate_guard.sh: unknown option $1" >&2; exit 2 ;;
        esac
    done
    if [[ -z "$mapper" || -z "$flag_dir" || -z "$sample" || -z "$min" || -z "$report" || $# -eq 0 ]]; then
        echo "mapping_rate_guard.sh run: missing required option or command" >&2
        exit 2
    fi
    label="${label:-$sample}"
    mkdir -p "$flag_dir"

    # fail REASON_STAGE READS PCT -- record the verdict (if first), explain, exit 42
    fail() {
        verdict_set "$flag_dir" "$sample" fail "$mapper" "$label" "$1" "$2" "$3" "$min" || true
        fail_banner "$mapper" "$label" "$1" "$2" "$3" "$min"
        exit "$EXIT_LOW_MAPPING"
    }

    # judge STAGE "READS PCT" -- decide on one measurement
    judge() {
        local stage="$1" reads pct
        read -r reads pct <<<"$2"
        echo "[mapping_rate_guard] ${label}: ${pct}% $(rate_name "$mapper") after ${reads} reads (${stage}, threshold ${min}%)"
        if pct_below "$pct" "$min"; then
            kill_mapper
            fail "$stage" "$reads" "$pct"
        fi
        verdict_set "$flag_dir" "$sample" pass "$mapper" "$label" "$stage" "$reads" "$pct" "$min" || true
    }

    if [[ "$(verdict_get "$flag_dir" "$sample")" == "fail" ]]; then
        fail other-mapper 0 NA
    fi

    "$@" &
    local pid=$!
    kill_mapper() { kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null; }
    trap 'kill "$pid" 2>/dev/null' EXIT
    trap 'kill "$pid" 2>/dev/null; exit 143' TERM INT

    # Monitor while the mapper runs. `sleep & wait` keeps the TERM trap responsive.
    while kill -0 "$pid" 2>/dev/null; do
        sleep "$poll" & wait $! 2>/dev/null
        kill -0 "$pid" 2>/dev/null || break

        case "$(verdict_get "$flag_dir" "$sample")" in
            fail) kill_mapper; fail other-mapper 0 NA ;;
            pass) continue ;;
        esac

        if [[ "$mapper" == "star" && -n "$progress" && "$check_reads" -gt 0 ]]; then
            local measure
            measure=$(star_progress "$progress")
            if [[ -n "$measure" ]] && awk -v r="${measure%% *}" -v n="$check_reads" 'BEGIN { exit !(r + 0 >= n + 0) }'; then
                judge mid-run "$measure"
            fi
        fi
    done

    local rc=0
    wait "$pid" || rc=$?
    trap - EXIT TERM INT
    if [[ $rc -ne 0 ]]; then
        # The other mapper may have cancelled us between polls: report that, not a crash.
        [[ "$(verdict_get "$flag_dir" "$sample")" == "fail" ]] && fail other-mapper 0 NA
        exit "$rc"
    fi

    # Final check on the finished mapper's own summary, unless a verdict already exists.
    case "$(verdict_get "$flag_dir" "$sample")" in
        fail) fail other-mapper 0 NA ;;
        pass) exit 0 ;;
    esac
    local measure
    if [[ "$mapper" == "star" ]]; then measure=$(star_final "$report"); else measure=$(salmon_meta "$report"); fi
    if [[ -z "$measure" ]]; then
        echo "[mapping_rate_guard] WARNING: no mapping rate found in ${report}; check skipped for ${label}" >&2
        exit 0
    fi
    kill_mapper() { :; }
    judge final "$measure"
    exit 0
}

# --------------------------------------------------------------------------

cmd="${1:-}"; shift || true
case "$cmd" in
    run)           run "$@" ;;
    star-progress) star_progress "$1" ;;
    star-final)    star_final "$1" ;;
    salmon)        salmon_meta "$1" ;;
    verdict-get)   verdict_get "$1" "$2" ;;
    verdict-set)   verdict_set "$@" ;;
    *) sed -n '2,36p' "$0" | sed 's/^# \{0,1\}//'; exit 2 ;;
esac
