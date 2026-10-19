#!/usr/bin/env bash
# description: Exercise lib/BcaResources.groovy, which sizes memory from input size
#
# This code decides how much memory a task asks the scheduler for. A mistake here
# does not throw -- it produces a plausible number that is too small, and the task
# is killed halfway through a run. The cases below therefore pin the behaviour that
# has to hold no matter what: an unmeasured process must get exactly its label
# value, a malformed entry must never throw, a process holding a large resident
# reference must never be handed less memory than that reference alone needs, and an
# input must still be measurable once Nextflow has wrapped it as a staged path.
#
# The same file also divides an allocation up inside a task, as STAR's
# --limitBAMsortRAM, which fails the same silent way.
#
# Groovy is run from the Nextflow fat jar that Nextflow itself caches, so no extra
# toolchain is required -- but the check skips cleanly when neither is present.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

JAVA_BIN="${BCA_JAVA:-}"
JAR="${BCA_NEXTFLOW_JAR:-}"

usage() {
    cat <<EOF
Usage: tests/run_tests.sh dynamic_resources [-- OPTIONS]
       tests/checks/dynamic_resources.sh [OPTIONS]

Unit-test lib/BcaResources.groovy, which turns an input size into a memory request
and an allocation into STAR's --limitBAMsortRAM.

Options:
  --java PATH     Java binary (default: \$BCA_JAVA, else java on PATH).
  --jar PATH      Nextflow fat jar providing Groovy and the config parser
                  (default: \$BCA_NEXTFLOW_JAR, else the newest under
                  \$NXF_HOME/framework/*/nextflow-*-one.jar).
  -h, --help      Show this message.

Cases:
  no_spec           an unmeasured process gets its label value, unchanged
  malformed_spec    partial, empty and undefined entries fall back, never throw
  scaling           the request tracks the input size between floor and cap
  cap_above_label   cap_gb can exceed the label, for a process near its tier limit
  resident_floor    a large resident reference raises the floor
  retry_ladder      task.attempt still multiplies the request
  staged_input      an input arriving as a staged TaskPath is measured, not read as 0
  staged_resident_floor  the resident floor survives that same trip
  bamsort_derived   STAR's BAM sort budget is the allocation minus index and slack
  bamsort_pinned    params.star_limitBAMsortRAM is used verbatim
  bamsort_no_memory no allocation to derive from leaves STAR its own default
  bamsort_unmeasured  an unmeasurable index falls back to a fraction, not a guess
  bamsort_floor     an index larger than the allocation floors at 1 GB
  config_shape      the dynamic_memory map in nextflow.config parses and is complete
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --java) JAVA_BIN="${2:?--java needs a value}"; shift 2 ;;
        --jar)  JAR="${2:?--jar needs a value}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

log_header "dynamic_resources"

[[ -n "$JAVA_BIN" ]] || { have_cmd java && JAVA_BIN="java"; }
if [[ -z "$JAR" ]]; then
    JAR="$(ls -1 "${NXF_HOME:-$HOME/.nextflow}"/framework/*/nextflow-*-one.jar 2>/dev/null | sort -V | tail -1 || true)"
fi

if [[ -z "$JAVA_BIN" ]] || ! "$JAVA_BIN" -version >/dev/null 2>&1; then
    record SKIP "java" "no java available"
    finish_check
    exit 0
fi
if [[ -z "$JAR" || ! -f "$JAR" ]]; then
    record SKIP "groovy" "no Nextflow jar found; run the pipeline once, or pass --jar"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/dynamic_resources"
mkdir -p "$WORKDIR/small_input" "$WORKDIR/big_reference"
# 256 MB of input, 2 GB of resident reference.
head -c 268435456 /dev/zero >"$WORKDIR/small_input/data.mtx" 2>/dev/null || true
head -c 2147483648 /dev/zero >"$WORKDIR/big_reference/SA" 2>/dev/null || true

cat >"$WORKDIR/Spec.groovy" <<'GROOVY'
import java.nio.file.Paths
import nextflow.config.ConfigParserFactory
import nextflow.file.FileHolder
import nextflow.processor.TaskPath
import nextflow.util.MemoryUnit

int failures = 0
void check(String name, boolean ok, String detail) {
    println((ok ? "PASS" : "FAIL") + "\t" + name + "\t" + detail)
}

def work    = Paths.get(System.getProperty("bca.work"))
def root    = System.getProperty("bca.root")
def input   = work.resolve("small_input")      // 256 MB
def bigRef  = work.resolve("big_reference")    // 2 GB
def spec    = [exponent: 0.89d, ref_gb: 5, mem_gb: 11, floor_gb: 1, cap_gb: 48]

// An unmeasured process must be indistinguishable from before the change.
check("no_spec", BcaResources.scaledMemory(null, input, 1, 12) == "12 GB",
      BcaResources.scaledMemory(null, input, 1, 12) + " (want 12 GB)")

// Config gives ConfigObject, whose missing keys are empty objects rather than
// null -- the guard has to survive that, and every malformed shape.
def parsed = ConfigParserFactory.create()
        .parse("params { dynamic_memory = [ PRESENT: [ exponent: 0.5, ref_gb: 1, mem_gb: 4 ], PARTIAL: [ exponent: 0.5 ] ] }")
        .params.dynamic_memory
def malformed = [
    BcaResources.scaledMemory(parsed?.ABSENT,  input, 1, 12),
    BcaResources.scaledMemory(parsed?.PARTIAL, input, 1, 12),
    BcaResources.scaledMemory([:],             input, 1, 12),
    BcaResources.scaledMemory(spec, Paths.get("does-not-exist"), 1, 12),
]
check("malformed_spec", malformed.every { it == "12 GB" }, malformed.join(", ") + " (all want 12 GB)")

// The request must track the input, and stay inside [floor, cap].
def small = BcaResources.scaledMemory(spec, input, 1, 12)          // 0.25 GB in
def wide  = [exponent: 0.89d, ref_gb: 0.25, mem_gb: 4, floor_gb: 1, cap_gb: 48]
def mid   = BcaResources.scaledMemory(wide, input, 1, 12)          // anchored at the input
check("scaling", small == "1 GB" && mid == "4 GB",
      "small=" + small + " (want 1 GB), anchored=" + mid + " (want 4 GB)")

// A process already near its tier limit needs room above it, or the cap simply
// reinstates the failure the scaling was meant to avoid.
def forcing = [exponent: 0.89d, ref_gb: 0.001, mem_gb: 40, floor_gb: 1, cap_gb: 48]
def above   = BcaResources.scaledMemory(forcing, input, 1, 12)
def capped  = BcaResources.scaledMemory([exponent: 0.89d, ref_gb: 0.001, mem_gb: 40, floor_gb: 1], input, 1, 12)
check("cap_above_label", above == "48 GB" && capped == "12 GB",
      "with cap_gb=" + above + " (want 48 GB), without=" + capped + " (want 12 GB)")

// A 2 GB resident reference must lift the floor above the scaled term.
def tiny      = [exponent: 0.89d, ref_gb: 5, mem_gb: 11, floor_gb: 1, cap_gb: 48]
def noRef     = BcaResources.scaledMemory(tiny, input, 1, 48)
def withRef   = BcaResources.scaledMemory(tiny, input, 1, 48, bigRef)
check("resident_floor", (withRef.replace(" GB","") as int) > (noRef.replace(" GB","") as int) &&
      withRef == "3 GB", "without=" + noRef + ", with 2 GB reference=" + withRef + " (want 3 GB)")

// The retry ladder has to survive, or a task killed for memory never escalates.
def a1 = BcaResources.scaledMemory(spec, input, 1, 12)
def a2 = BcaResources.scaledMemory(spec, input, 2, 12)
check("retry_ladder", (a2.replace(" GB","") as int) == 2 * (a1.replace(" GB","") as int),
      "attempt1=" + a1 + " attempt2=" + a2)

// Inside a process, Nextflow does not hand over the path it read the input from:
// it hands over a TaskPath, whose provider is not the one that owns the file, so
// every java.nio.file.Files call on it throws. Measured through it, an input reads
// as 0 bytes and every request silently collapses onto its label -- which is the
// whole feature failing quietly. This is what the input actually looks like.
def staged      = new TaskPath(new FileHolder(bigRef))
def stagedBytes = BcaResources.totalBytes(staged)
check("staged_input", stagedBytes == BcaResources.totalBytes(bigRef) && stagedBytes == 2147483648L,
      "staged=" + stagedBytes + " (want 2147483648, the same as the unstaged path)")

// ... and the resident floor has to survive that same trip, since the genome index
// reaches STARSOLO_ALIGN as a staged directory and nothing else.
check("staged_resident_floor", BcaResources.scaledMemory(tiny, input, 1, 48, staged) == "3 GB",
      BcaResources.scaledMemory(tiny, input, 1, 48, staged) + " (want 3 GB)")

// --limitBAMsortRAM: the allocation, less the genome index STAR keeps resident and
// its own working memory. A value too large here is not caught by the scheduler --
// STAR simply allocates past the cgroup limit and the task is killed.
def GB       = 1024L * 1024 * 1024
def alloc    = new MemoryUnit("32 GB")
def derived  = BcaResources.bamSortRam(0, alloc, bigRef)          // 32 - 2 - 4 slack
check("bamsort_derived", derived.bytes == 26L * GB,
      derived.bytes + " (want " + (26L * GB) + "); note: " + derived.note)

// A pinned param wins outright: it is how a user answers a STAR error that names
// the exact figure it wants.
check("bamsort_pinned", BcaResources.bamSortRam(123456789L, alloc, bigRef).bytes == 123456789L,
      BcaResources.bamSortRam(123456789L, alloc, bigRef).note)

// Nothing to derive from must leave STAR its own behaviour, not a guess.
check("bamsort_no_memory", BcaResources.bamSortRam(0, null, bigRef).bytes == 0L,
      BcaResources.bamSortRam(0, null, bigRef).note)

// An index that cannot be measured must not be assumed small.
def unknown = BcaResources.bamSortRam(0, alloc, Paths.get("does-not-exist"))
check("bamsort_unmeasured", unknown.bytes == (32L * GB * 60L).intdiv(100L),
      unknown.bytes + " (want 60% of 32 GB); note: " + unknown.note)

// An index too big for the allocation leaves nothing to sort in: ask for a modest
// buffer rather than a negative one and let STAR report the real shortfall.
def squeezed = BcaResources.bamSortRam(0, new MemoryUnit("2 GB"), bigRef)
check("bamsort_floor", squeezed.bytes == GB, squeezed.bytes + " (want " + GB + "); note: " + squeezed.note)

// Every entry shipped in nextflow.config must be complete, or it silently does
// nothing at all.
def text = new File(root, "nextflow.config").text
def start = text.indexOf("dynamic_memory = [")
def depth = 0, end = start
while (end < text.length()) {
    if (text.charAt(end) == ('[' as char)) depth++
    else if (text.charAt(end) == (']' as char)) { depth--; if (depth == 0) break }
    end++
}
def shipped = ConfigParserFactory.create()
        .parse("params {\n" + text.substring(start, end + 1) + "\n}").params.dynamic_memory
def incomplete = shipped.findAll { k, v -> !(v.exponent != null && v.ref_gb && v.mem_gb) }.keySet()
check("config_shape", shipped.size() > 0 && incomplete.isEmpty(),
      shipped.size() + " entries: " + shipped.keySet().join(", ") +
      (incomplete ? "; INCOMPLETE: " + incomplete.join(", ") : ""))
GROOVY

# The JVM needs native paths and the platform's classpath separator, which under
# Git Bash on Windows is neither the shell's /c/... form nor a colon.
to_native() {
    if have_cmd cygpath; then cygpath -w "$1"; else printf '%s' "$1"; fi
}
CP_SEP=":"
case "$(uname -s 2>/dev/null)" in
    CYGWIN*|MINGW*|MSYS*) CP_SEP=";" ;;
esac
CLASSPATH_ARG="$(to_native "$JAR")${CP_SEP}$(to_native "$PROJECT_ROOT/lib")"

OUT="$WORKDIR/results.tsv"
if "$JAVA_BIN" -cp "$CLASSPATH_ARG" \
        -Dbca.work="$(to_native "$WORKDIR")" -Dbca.root="$(to_native "$PROJECT_ROOT")" \
        groovy.ui.GroovyMain "$(to_native "$WORKDIR/Spec.groovy")" 2>"$WORKDIR/stderr.log" \
        | grep -E '^(PASS|FAIL)' >"$OUT"; then
    while IFS=$'\t' read -r status name detail; do
        record "$status" "$name" "$detail"
    done <"$OUT"
    [[ -s "$OUT" ]] || record FAIL "groovy" "no results produced"
else
    record FAIL "groovy" "$(tail -3 "$WORKDIR/stderr.log" 2>/dev/null | tr '\n' ' ')"
fi

finish_check
