#!/usr/bin/env bash
# description: Check the mapping-rate guard that cancels low-mapping samples
#
# STARSOLO_ALIGN and ALEVIN_FRY run their mapper through bin/mapping_rate_guard.sh,
# which reads the mapping rate from STAR's Log.progress.out (mid-run) and
# Log.final.out, or from salmon's meta_info.json, and exits with code 42 when a
# sample maps below the threshold. Both mappers of a sample share a verdict file,
# and the first verdict wins.
#
# The mappers are replaced by small shell scripts that write synthetic logs, so
# no STAR, salmon, sequencing data or cluster is needed.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

GUARD="$PROJECT_ROOT/bin/mapping_rate_guard.sh"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh mapping_rate_guard [-- OPTIONS]
       tests/checks/mapping_rate_guard.sh [OPTIONS]

Validate bin/mapping_rate_guard.sh against synthetic mapper logs.

Options:
  --keep          Keep the generated fixtures for inspection.
  -h, --help      Show this message.

Cases:
  parse_progress     last data row of Log.progress.out is read
  parse_progress_hdr a Log.progress.out with only headers yields nothing
  parse_final        reads and unique % from Log.final.out
  parse_salmon       num_processed and percent_mapped from meta_info.json
  verdict_first      the first verdict wins, a second one is refused
  midrun_kill        a running STAR below the threshold is killed, exit 42
  midrun_pass        a passing mid-run check lets STAR finish, exit 0
  final_fail         a low Log.final.out after STAR finishes, exit 42
  salmon_fail        a low percent_mapped after salmon finishes, exit 42
  other_mapper       a recorded fail cancels the mapper before it starts
  other_mapper_live  a fail recorded mid-run by the other mapper kills this one
  pass_skips_check   a recorded pass skips this mapper's own check
  mapper_error       a mapper crash keeps its own exit code
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep)    KEEP=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

log_header "mapping_rate_guard"

WORKDIR="$BCA_TEST_LOGDIR/mapping_rate_guard"
rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

progress_header() {
    cat <<'EOF'
           Time    Speed        Read     Read   Mapped   Mapped   Mapped   Mapped Unmapped Unmapped Unmapped Unmapped
                    M/hr      number   length   unique   length   MMrate    multi   multi+       MM    short    other
EOF
}

# progress_row READS PCT
progress_row() {
    printf 'Sep 30 10:01:23     200.3 %11s      132 %7s%%    130.2     0.4%%     5.2%%     0.1%%     0.0%%     9.2%%     0.4%%\n' "$1" "$2"
}

# final_log FILE READS PCT
final_log() {
    cat >"$1" <<EOF
                                 Started job on |	Sep 30 10:00:00
                          Number of input reads |	$2
                      Average input read length |	132
                                    UNIQUE READS:
                   Uniquely mapped reads number |	123
                        Uniquely mapped reads % |	$3%
                          Average mapped length |	130.20
EOF
}

# salmon_meta FILE READS PCT
salmon_meta() {
    mkdir -p "$(dirname "$1")"
    cat >"$1" <<EOF
{
    "salmon_version": "1.10.1",
    "num_processed": $2,
    "num_mapped": 1000,
    "percent_mapped": $3,
    "library_types": ["ISR"]
}
EOF
}

# Fake STAR: writes progress rows once per second, then optionally Log.final.out.
#   fake_star DIR PREFIX PCT ROWS FINAL_PCT
cat >"$WORKDIR/fake_star.sh" <<'EOF'
#!/usr/bin/env bash
dir="$1"; prefix="$2"; pct="$3"; rows="$4"; final_pct="$5"
cd "$dir"
{
cat <<'H'
           Time    Speed        Read     Read   Mapped   Mapped   Mapped   Mapped Unmapped Unmapped Unmapped Unmapped
                    M/hr      number   length   unique   length   MMrate    multi   multi+       MM    short    other
H
} >"${prefix}Log.progress.out"
for i in $(seq 1 "$rows"); do
    printf 'Sep 30 10:0%d:00     200.3 %11s      132 %7s%%    130.2     0.4%%     5.2%%     0.1%%     0.0%%     9.2%%     0.4%%\n' \
        "$((i % 10))" "$((i * 1000000))" "$pct" >>"${prefix}Log.progress.out"
    sleep 1
done
touch "${prefix}finished"
[[ -n "$final_pct" ]] && cat >"${prefix}Log.final.out" <<F
                          Number of input reads |	$((rows * 1000000))
                        Uniquely mapped reads % |	${final_pct}%
F
exit 0
EOF
chmod +x "$WORKDIR/fake_star.sh" "$GUARD" 2>/dev/null || true

# run_guard CASE_DIR MAPPER [extra guard options...] -- COMMAND...
# Sets RC and captures stdout/stderr in CASE_DIR/out.log.
run_guard() {
    local dir="$1" mapper="$2"; shift 2
    RC=0
    ( cd "$dir" && bash "$GUARD" run --mapper "$mapper" --flag-dir "$dir/flags" \
        --sample S1 --label "S1_${mapper}" --min 10 --poll 1 "$@" ) >"$dir/out.log" 2>&1 || RC=$?
}

case_dir() { local d="$WORKDIR/$1"; mkdir -p "$d"; echo "$d"; }

expect() {  # expect NAME CONDITION_RESULT DETAIL
    if [[ "$2" -eq 0 ]]; then record PASS "$1"; else record FAIL "$1" "$3"; fi
}

# --------------------------------------------------------------------------
# Parsers
# --------------------------------------------------------------------------

log_step "parsers"

f="$WORKDIR/progress.out"
{ progress_header; progress_row 1000000 85.0; progress_row 25000000 4.3; echo "ALL DONE!"; } >"$f"
got=$(bash "$GUARD" star-progress "$f")
[[ "$got" == "25000000 4.3" ]] && r=0 || r=1; expect parse_progress $r "got '$got'"

progress_header >"$f"
got=$(bash "$GUARD" star-progress "$f")
[[ -z "$got" ]] && r=0 || r=1; expect parse_progress_hdr $r "got '$got'"

final_log "$WORKDIR/Log.final.out" 5000000 61.25
got=$(bash "$GUARD" star-final "$WORKDIR/Log.final.out")
[[ "$got" == "5000000 61.25" ]] && r=0 || r=1; expect parse_final $r "got '$got'"

salmon_meta "$WORKDIR/meta_info.json" 7000000 3.5
got=$(bash "$GUARD" salmon "$WORKDIR/meta_info.json")
[[ "$got" == "7000000 3.5" ]] && r=0 || r=1; expect parse_salmon $r "got '$got'"

d=$(case_dir verdict_first)
ok=0
bash "$GUARD" verdict-set "$d" S1 fail star S1_starsolo mid-run 20000000 4.3 10 || ok=1
bash "$GUARD" verdict-set "$d" S1 pass salmon S1_alevinfry final 20000000 50 10 2>/dev/null && ok=1
[[ "$(bash "$GUARD" verdict-get "$d" S1)" == "fail" && "$(bash "$GUARD" verdict-get "$d" S2)" == "none" ]] || ok=1
expect verdict_first $ok "$(cat "$d"/*.verdict 2>/dev/null)"

# --------------------------------------------------------------------------
# Wrapped mapper runs
# --------------------------------------------------------------------------

log_step "wrapped runs"

# STAR at 4% for 60 rows (60 s); the check at 2M reads must kill it within seconds.
d=$(case_dir midrun_kill)
start=$SECONDS
run_guard "$d" star --check-reads 2000000 --progress "$d/X_Log.progress.out" --report "$d/X_Log.final.out" \
    -- "$WORKDIR/fake_star.sh" "$d" X_ 4.3 60 4.3
elapsed=$((SECONDS - start))
ok=0
[[ $RC -eq 42 ]] || ok=1
[[ $elapsed -lt 20 ]] || ok=1
[[ ! -e "$d/X_finished" ]] || ok=1
grep -q "MAPPING RATE CHECK FAILED" "$d/out.log" || ok=1
grep -q "Kraken2" "$d/out.log" || ok=1
[[ "$(cut -f2,5 "$d/flags/S1.verdict")" == $'fail\tmid-run' ]] || ok=1
expect midrun_kill $ok "rc=$RC elapsed=${elapsed}s"

d=$(case_dir midrun_pass)
run_guard "$d" star --check-reads 2000000 --progress "$d/X_Log.progress.out" --report "$d/X_Log.final.out" \
    -- "$WORKDIR/fake_star.sh" "$d" X_ 85.0 4 85.0
ok=0
[[ $RC -eq 0 && -e "$d/X_finished" ]] || ok=1
[[ "$(cut -f2,5 "$d/flags/S1.verdict")" == $'pass\tmid-run' ]] || ok=1
expect midrun_pass $ok "rc=$RC verdict=$(cat "$d/flags/S1.verdict" 2>/dev/null)"

# Too few reads for the mid-run check: decided on Log.final.out.
d=$(case_dir final_fail)
run_guard "$d" star --check-reads 20000000 --progress "$d/X_Log.progress.out" --report "$d/X_Log.final.out" \
    -- "$WORKDIR/fake_star.sh" "$d" X_ 85.0 2 4.0
ok=0
[[ $RC -eq 42 ]] || ok=1
[[ "$(cut -f2,5 "$d/flags/S1.verdict")" == $'fail\tfinal' ]] || ok=1
expect final_fail $ok "rc=$RC"

d=$(case_dir salmon_fail)
run_guard "$d" salmon --report "$d/run/aux_info/meta_info.json" \
    -- bash -c "sleep 1; mkdir -p '$d/run/aux_info'; cat '$WORKDIR/meta_info.json' > '$d/run/aux_info/meta_info.json'"
ok=0
[[ $RC -eq 42 ]] || ok=1
grep -q "alevin-fry" "$d/out.log" || ok=1
expect salmon_fail $ok "rc=$RC"

d=$(case_dir other_mapper)
bash "$GUARD" verdict-set "$d/flags" S1 fail star S1_starsolo mid-run 20000000 4.3 10
run_guard "$d" salmon --report "$d/meta_info.json" -- touch "$d/started"
ok=0
[[ $RC -eq 42 && ! -e "$d/started" ]] || ok=1
grep -q "other mapper" "$d/out.log" || ok=1
expect other_mapper $ok "rc=$RC"

d=$(case_dir other_mapper_live)
( sleep 2; bash "$GUARD" verdict-set "$d/flags" S1 fail star S1_starsolo mid-run 20000000 4.3 10 ) &
start=$SECONDS
run_guard "$d" salmon --report "$d/meta_info.json" -- bash -c "sleep 60; touch '$d/finished'"
elapsed=$((SECONDS - start))
wait
ok=0
[[ $RC -eq 42 && $elapsed -lt 20 && ! -e "$d/finished" ]] || ok=1
expect other_mapper_live $ok "rc=$RC elapsed=${elapsed}s"

d=$(case_dir pass_skips_check)
bash "$GUARD" verdict-set "$d/flags" S1 pass star S1_starsolo final 20000000 60 10
salmon_meta "$d/meta_info.json" 7000000 3.5
run_guard "$d" salmon --report "$d/meta_info.json" -- true
[[ $RC -eq 0 ]] && r=0 || r=1; expect pass_skips_check $r "rc=$RC"

d=$(case_dir mapper_error)
run_guard "$d" star --check-reads 0 --progress "$d/none" --report "$d/none" -- bash -c "exit 3"
[[ $RC -eq 3 ]] && r=0 || r=1; expect mapper_error $r "rc=$RC"

finish_check
