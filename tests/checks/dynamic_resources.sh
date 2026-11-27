#!/usr/bin/env bash
# description: Exercise lib/BcaResources.groovy, which sizes memory, time and cpus from input size
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
  nested_inputs     a list of multi-file inputs is measured in full, not read as 0
  files_ending_with only a directory's matching files are measured (Kraken2 *.k2d)
  files_in          only the named files of a staged directory are measured (STARsolo raw/)
  time_no_spec      an unmeasured process gets exactly its time limit
  time_malformed    partial, empty and undefined time entries fall back, never throw
  time_scaling      the time limit tracks the input size between floor and cap
  time_floor_default without floor_h the label is the floor, so a limit never shrinks
  time_cap          cap_h bounds the first attempt; without it the label does
  time_retry_ladder task.attempt still multiplies the time limit
  time_parses       every returned limit parses as a Nextflow Duration
  time_staged_input a staged TaskPath input is measured for time too
  cpus_off          tieredCpus gives the base count while its threshold is unset
  cpus_tiered       tieredCpus steps up above the threshold, times task.attempt
  bamsort_derived   STAR's BAM sort budget is the allocation minus index and slack
  bamsort_pinned    params.star_limitBAMsortRAM is used verbatim
  bamsort_no_memory no allocation to derive from leaves STAR its own default
  bamsort_unmeasured  an unmeasurable index falls back to a fraction, not a guess
  bamsort_floor     an index larger than the allocation floors at 1 GB
  retry_by_exit     a retry doubles time after exit 140, memory after 137, both otherwise
  retry_via_task    attemptMemory/attemptTime, the modules' form, read attempt and trace from the task
  as_strings        "8 GB" and "131m" read back as bytes and milliseconds
  prior_time        a GeneExt re-mapping's limit is the earlier pass's time x prior_factor
  star_elapsed      STAR's run time is read from Log.final.out, also across New Year
  config_shape      the dynamic_memory, dynamic_time and dynamic_cpus maps in
                    nextflow.config parse and are complete
  star_memory_keys  STARSOLO_ALIGN and STARSOLO_ALIGN_BAM_ONLY are both shipped
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
mkdir -p "$WORKDIR/small_input" "$WORKDIR/big_reference" "$WORKDIR/kraken_db/library"
# 256 MB of input, 2 GB of resident reference.
head -c 268435456 /dev/zero >"$WORKDIR/small_input/data.mtx" 2>/dev/null || true
head -c 2147483648 /dev/zero >"$WORKDIR/big_reference/SA" 2>/dev/null || true
# A Kraken2 database built in place: 3 MB of *.k2d, plus build leftovers that are not loaded.
head -c 2097152 /dev/zero >"$WORKDIR/kraken_db/hash.k2d"
head -c 1048576 /dev/zero >"$WORKDIR/kraken_db/taxo.k2d"
head -c 4194304 /dev/zero >"$WORKDIR/kraken_db/seqid2taxid.map"
head -c 8388608 /dev/zero >"$WORKDIR/kraken_db/library/library.fna"
# A STARsolo Solo.out: 3 MB in the files a step reads, a 4 MB matrix next to them it does not.
mkdir -p "$WORKDIR/solo_out/Gene/raw"
head -c 2097152 /dev/zero >"$WORKDIR/solo_out/Gene/raw/matrix.mtx"
head -c 1048576 /dev/zero >"$WORKDIR/solo_out/Gene/raw/barcodes.tsv"
head -c 4194304 /dev/zero >"$WORKDIR/solo_out/Gene/raw/UniqueAndMult-EM.mtx"
# STAR Log.final.out heads: a 3h14m run (hoihon01_H1H13_1, genoa64), one across New Year,
# and one cut off before "Finished on".
mkdir -p "$WORKDIR/star_logs"
printf '%s\n' \
    $'                                 Started job on |\tNov 23 12:08:36' \
    $'                             Started mapping on |\tNov 23 12:09:02' \
    $'                                    Finished on |\tNov 23 15:22:41' \
    $'       Mapping speed, Million of reads per hour |\t43.21' \
    $'                          Number of input reads |\t139000000' \
    >"$WORKDIR/star_logs/H1H13_Log.final.out"
printf '%s\n' \
    $'                                 Started job on |\tDec 31 23:30:00' \
    $'                                    Finished on |\tJan 01 00:45:30' \
    >"$WORKDIR/star_logs/newyear_Log.final.out"
printf '%s\n' \
    $'                                 Started job on |\tNov 23 12:08:36' \
    >"$WORKDIR/star_logs/cut_Log.final.out"

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

// A path input that received a collected list reaches the directive as a list of its
// own (a BlankSeparatedList), so GENERATE_DASHBOARD's inputs are a list of lists.
def multi       = new nextflow.util.BlankSeparatedList([staged, input.resolve("data.mtx")])
def nestedBytes = BcaResources.totalBytes([multi, [input]])
check("nested_inputs", nestedBytes == 2147483648L + 2L * 268435456L,
      "nested=" + nestedBytes + " (want " + (2147483648L + 2L * 268435456L) + ")")

// Only the *.k2d files of a Kraken2 database are loaded; the library/ a database built
// in place keeps, and its maps, must not count. The directory arrives staged.
def k2d      = BcaResources.filesEndingWith(new TaskPath(new FileHolder(work.resolve("kraken_db"))), ".k2d")
def k2dBytes = BcaResources.totalBytes(k2d)
def noDir    = BcaResources.filesEndingWith(Paths.get("does-not-exist"), ".k2d")
check("files_ending_with", k2d.size() == 2 && k2dBytes == 3145728L && noDir.isEmpty(),
      k2d.size() + " files, " + k2dBytes + " bytes (want 2 files, 3145728 bytes); missing dir gives " + noDir.size())

// Only the named files of a staged directory are measured: STARsolo's raw/ holds every
// multi-gene matrix next to matrix.mtx, and only the analysed one is loaded. Absent
// names, empty names and subdirectories are left out; nested paths resolve.
def soloDir  = new TaskPath(new FileHolder(work.resolve("solo_out")))
def named    = BcaResources.filesIn(soloDir,
        ["Gene/raw/matrix.mtx", "Gene/raw/barcodes.tsv", "Velocyto/raw/spliced.mtx", "", null, "Gene/raw"])
def namedBytes = BcaResources.totalBytes(named)
def noSolo   = BcaResources.filesIn(Paths.get("does-not-exist"), ["matrix.mtx"])
check("files_in", named.size() == 2 && namedBytes == 3145728L && noSolo.isEmpty(),
      named.size() + " files, " + namedBytes + " bytes (want 2 files, 3145728 bytes; the 4 MB " +
      "UniqueAndMult-EM.mtx left out); missing dir gives " + noSolo.size())

// Time: an unmeasured process keeps exactly its limit, in a form Nextflow parses.
check("time_no_spec", BcaResources.scaledTime(null, input, 1, 0.75) == "45m" &&
      BcaResources.scaledTime(null, input, 1, 7) == "420m",
      BcaResources.scaledTime(null, input, 1, 0.75) + ", " + BcaResources.scaledTime(null, input, 1, 7) +
      " (want 45m, 420m)")

def parsedTime = ConfigParserFactory.create()
        .parse("params { dynamic_time = [ PARTIAL: [ exponent: 1.0 ] ] }")
        .params.dynamic_time
def tspec = [exponent: 1.0d, ref_gb: 0.25, time_h: 1, floor_h: 0.5, cap_h: 24]
def tmalformed = [
    BcaResources.scaledTime(parsedTime?.ABSENT,  input, 1, 2),
    BcaResources.scaledTime(parsedTime?.PARTIAL, input, 1, 2),
    BcaResources.scaledTime([:],                 input, 1, 2),
    BcaResources.scaledTime(tspec, Paths.get("does-not-exist"), 1, 2),
]
check("time_malformed", tmalformed.every { it == "120m" }, tmalformed.join(", ") + " (all want 120m)")

// The limit tracks the input: anchored at the input it is time_h, at 8x the anchor 8x.
def atRef   = BcaResources.scaledTime(tspec, input, 1, 0.5)
def bigger  = BcaResources.scaledTime([exponent: 1.0d, ref_gb: 0.03125, time_h: 1, floor_h: 0.5, cap_h: 24], input, 1, 0.5)
check("time_scaling", atRef == "60m" && bigger == "480m",
      "anchored=" + atRef + " (want 60m), 8x the anchor=" + bigger + " (want 480m)")

// Without floor_h the label is the floor: an entry can lengthen a limit, never shorten it.
def shortRun = [exponent: 1.0d, ref_gb: 1, time_h: 0.1, cap_h: 24]
def explicit = [exponent: 1.0d, ref_gb: 1, time_h: 0.1, floor_h: 0.25, cap_h: 24]
check("time_floor_default", BcaResources.scaledTime(shortRun, input, 1, 2) == "120m" &&
      BcaResources.scaledTime(explicit, input, 1, 2) == "15m",
      "no floor_h=" + BcaResources.scaledTime(shortRun, input, 1, 2) + " (want 120m), floor_h 0.25=" +
      BcaResources.scaledTime(explicit, input, 1, 2) + " (want 15m)")

// cap_h bounds the first attempt; without it the label does.
def huge = [exponent: 1.0d, ref_gb: 0.001, time_h: 10, floor_h: 1]
check("time_cap", BcaResources.scaledTime(huge + [cap_h: 48], input, 1, 2) == "2880m" &&
      BcaResources.scaledTime(huge, input, 1, 2) == "120m",
      "cap_h 48=" + BcaResources.scaledTime(huge + [cap_h: 48], input, 1, 2) + " (want 2880m), none=" +
      BcaResources.scaledTime(huge, input, 1, 2) + " (want 120m)")

check("time_retry_ladder", BcaResources.scaledTime(tspec, input, 2, 0.5) == "120m" &&
      BcaResources.scaledTime(null, input, 2, 0.75) == "90m",
      "attempt 2: " + BcaResources.scaledTime(tspec, input, 2, 0.5) + " (want 120m), unmeasured " +
      BcaResources.scaledTime(null, input, 2, 0.75) + " (want 90m)")

// Every value has to be a Duration, or the directive fails when the task is submitted.
def durations = ["45m", "420m", "2880m", BcaResources.scaledTime(tspec, input, 1, 1.0d / 3)]
def parsedDurations = durations.collect { new nextflow.util.Duration(it).toMinutes() }
check("time_parses", parsedDurations == [45L, 420L, 2880L, 60L],
      durations.join(", ") + " -> " + parsedDurations.join(", ") + " minutes")

check("time_staged_input",
      BcaResources.scaledTime([exponent: 1.0d, ref_gb: 2, time_h: 1, floor_h: 0.1, cap_h: 24], staged, 1, 1) == "60m",
      BcaResources.scaledTime([exponent: 1.0d, ref_gb: 2, time_h: 1, floor_h: 0.1, cap_h: 24], staged, 1, 1) +
      " (want 60m: the 2 GB staged reference at its own anchor)")

// CPUs: off while above_gb is unset (including absent in a ConfigObject), on above it.
def parsedCpus = ConfigParserFactory.create()
        .parse("params { dynamic_cpus { NOTHRESHOLD { cpus = 16 } } }").params.dynamic_cpus
def off = [
    BcaResources.tieredCpus(null, input, 1, 8),
    BcaResources.tieredCpus([above_gb: null, cpus: 16], input, 1, 8),
    BcaResources.tieredCpus(parsedCpus?.NOTHRESHOLD, input, 1, 8),
    BcaResources.tieredCpus([above_gb: 0.1, cpus: 16], Paths.get("does-not-exist"), 1, 8),
]
check("cpus_off", off == [8, 8, 8, 8], off.join(", ") + " (all want 8)")
def cpuBelow = BcaResources.tieredCpus([above_gb: 1, cpus: 16], input, 1, 8)
def cpuAbove = BcaResources.tieredCpus([above_gb: 0.1, cpus: 16], input, 1, 8)
def cpuRetry = BcaResources.tieredCpus([above_gb: 0.1, cpus: 16], input, 2, 8)
check("cpus_tiered", cpuBelow == 8 && cpuAbove == 16 && cpuRetry == 32,
      "0.25 GB under 1 GB=" + cpuBelow + " (want 8), over 0.1 GB=" + cpuAbove + " (want 16), attempt 2=" +
      cpuRetry + " (want 32)")

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

// Retries grow only what the failure called for, from the previous attempt's request.
// The previous attempt arrives as task.previousTrace (a TraceRecord: exit Integer,
// memory bytes, time milliseconds); a Map stands in for it here.
def GiB = 1024L * 1024 * 1024
def prevTrace = { int exit, long memGb, long minutes -> [exit: exit, memory: memGb * GiB, time: minutes * 60000L] }
def memAfter = { prev, int attempt -> new MemoryUnit(BcaResources.retryMemory(prev, attempt, "11 GB")).toGiga() }
def timeAfter = { prev, int attempt -> new nextflow.util.Duration(BcaResources.retryTime(prev, attempt, "131m")).toMinutes() }
def retries = [
    first:    [memAfter(null, 1), timeAfter(null, 1)],
    timeout:  [memAfter(prevTrace(140, 11, 131), 2), timeAfter(prevTrace(140, 11, 131), 2)],
    oom:      [memAfter(prevTrace(137, 11, 131), 2), timeAfter(prevTrace(137, 11, 131), 2)],
    other:    [memAfter(prevTrace(143, 11, 131), 2), timeAfter(prevTrace(143, 11, 131), 2)],
    // attempt 3 after an OOM (22 GB, 131m) then a timeout: memory stays doubled
    chained:  [memAfter(prevTrace(140, 22, 131), 3), timeAfter(prevTrace(140, 22, 131), 3)],
    // a previous trace without requests falls back to the plain ladder
    noTrace:  [memAfter([exit: 140], 2), timeAfter([exit: 140], 2)],
]
def wantRetries = [first: [11L, 131L], timeout: [11L, 262L], oom: [22L, 131L], other: [22L, 262L],
                   chained: [22L, 262L], noTrace: [22L, 262L]]
check("retry_by_exit", retries == wantRetries,
      retries.collect { k, v -> k + "=" + v[0] + "GB/" + v[1] + "m" }.join(", ") +
      " (want " + wantRetries.collect { k, v -> k + "=" + v[0] + "GB/" + v[1] + "m" }.join(", ") + ")")

// The module-facing form takes the task itself (a Map standing in for TaskConfig).
def viaTask = [
    BcaResources.attemptMemory([attempt: 1], "2 GB"),
    BcaResources.attemptMemory([attempt: 2, previousTrace: prevTrace(140, 2, 30)], "2 GB"),
    BcaResources.attemptMemory([attempt: 2, previousTrace: prevTrace(137, 2, 30)], "2 GB"),
    BcaResources.attemptTime([attempt: 1], "4h"),
    BcaResources.attemptTime([attempt: 2, previousTrace: prevTrace(140, 2, 240)], "4h"),
    BcaResources.attemptTime([attempt: 2, previousTrace: prevTrace(137, 2, 240)], "4h"),
]
def wantViaTask = ["2048 MB", "2048 MB", "4096 MB", "240m", "480m", "240m"]
check("retry_via_task", viaTask == wantViaTask, viaTask.join(", ") + " (want " + wantViaTask.join(", ") + ")")

// The strings this class writes read back as the values it meant.
def asStrings = [BcaResources.asBytes("8 GB"), BcaResources.asBytes("512 MB"), BcaResources.asBytes("1024"),
                 BcaResources.asBytes("auto"), BcaResources.asMillis("131m"), BcaResources.asMillis("4h"),
                 BcaResources.asMillis(new nextflow.util.Duration("2h")), BcaResources.asMillis("soon")]
def wantStrings = [8L * GiB, 512L * 1024 * 1024, 1024L, 0L, 131L * 60000, 4L * 3600000, 2L * 3600000, 0L]
check("as_strings", asStrings == wantStrings, asStrings.join(", ") + " (want " + wantStrings.join(", ") + ")")

// A GeneExt re-mapping's time from the earlier pass: minutes x prior_factor, within floor and cap.
def pspec = [exponent: 0.9d, ref_gb: 34.4, time_h: 4, floor_h: 1, cap_h: 24, prior_factor: 1.75d]
def parsedPrior = ConfigParserFactory.create()
        .parse("params { dynamic_time = [ NOPRIOR: [ exponent: 0.9, ref_gb: 34.4, time_h: 4 ] ] }").params.dynamic_time
def priors = [
    BcaResources.priorTime(pspec, 194, "131m"),                     // 194 x 1.75 = 339.5 -> 340m
    BcaResources.priorTime(pspec, 12, "60m"),                       // 21 min, floored at 1 h
    BcaResources.priorTime(pspec, 1000, "131m"),                    // 29 h, capped at 24 h
    BcaResources.priorTime(pspec, null, "131m"),                    // no earlier pass
    BcaResources.priorTime(pspec, 0, "131m"),
    BcaResources.priorTime(pspec - [prior_factor: 1.75d], 194, "131m"),
    BcaResources.priorTime(parsedPrior?.NOPRIOR, 194, "131m"),
    BcaResources.priorTime(null, 194, "131m"),
]
def wantPriors = ["340m", "60m", "1440m", "131m", "131m", "131m", "131m", "131m"]
check("prior_time", priors == wantPriors, priors.join(", ") + " (want " + wantPriors.join(", ") + ")")

// STAR's run time read from Log.final.out
def logs = Paths.get(System.getProperty("bca.work")).resolve("star_logs")
def elapsed = [
    BcaResources.starElapsedMinutes(logs.resolve("H1H13_Log.final.out")),               // 3h14m05s -> 195
    BcaResources.starElapsedMinutes(new TaskPath(new FileHolder(logs.resolve("H1H13_Log.final.out")))),
    BcaResources.starElapsedMinutes(logs.resolve("newyear_Log.final.out")),             // 1h15m30s -> 76
    BcaResources.starElapsedMinutes(logs.resolve("cut_Log.final.out")),
    BcaResources.starElapsedMinutes(logs.resolve("does-not-exist")),
]
check("star_elapsed", elapsed == [195, 195, 76, null, null],
      elapsed.join(", ") + " (want 195, 195, 76, null, null)")

// Every entry shipped in nextflow.config must be complete, or it silently does
// nothing at all.
def text = new File(root, "nextflow.config").text
def shippedMap = { String name ->
    def start = text.indexOf(name + " = [")
    if (start < 0) return [:]
    def depth = 0, end = start
    while (end < text.length()) {
        if (text.charAt(end) == ('[' as char)) depth++
        else if (text.charAt(end) == (']' as char)) { depth--; if (depth == 0) break }
        end++
    }
    ConfigParserFactory.create().parse("params {\n" + text.substring(start, end + 1) + "\n}").params[name]
}
def shipped = shippedMap("dynamic_memory")
def incomplete = shipped.findAll { k, v -> !(v.exponent != null && v.ref_gb && v.mem_gb) }.keySet()
def shippedTime = shippedMap("dynamic_time")
def incompleteTime = shippedTime.findAll { k, v -> !(v.exponent != null && v.ref_gb && v.time_h) }.keySet()
def shippedCpus = shippedMap("dynamic_cpus")
def badCpus = shippedCpus.findAll { k, v -> !(v.containsKey('above_gb') && v.cpus) }.keySet()
check("config_shape", shipped.size() > 0 && incomplete.isEmpty() && shippedTime.size() > 0 &&
      incompleteTime.isEmpty() && badCpus.isEmpty(),
      "dynamic_memory " + shipped.size() + ", dynamic_time " + shippedTime.size() + " (" +
      shippedTime.keySet().join(", ") + "), dynamic_cpus " + shippedCpus.size() + " entries" +
      (incomplete ? "; INCOMPLETE memory: " + incomplete.join(", ") : "") +
      (incompleteTime ? "; INCOMPLETE time: " + incompleteTime.join(", ") : "") +
      (badCpus ? "; INCOMPLETE cpus: " + badCpus.join(", ") : ""))

// The module picks STARSOLO_ALIGN_BAM_ONLY for a bam_only pass; a missing entry would
// silently fall back to the 32 GB label. The full pass must not ask for less.
def starFull = shipped.STARSOLO_ALIGN, starBam = shipped.STARSOLO_ALIGN_BAM_ONLY
check("star_memory_keys", starFull?.mem_gb && starBam?.mem_gb && starFull.floor_gb >= starBam.floor_gb &&
      starFull.mem_gb >= starBam.mem_gb && shippedTime.STARSOLO_ALIGN?.prior_factor instanceof Number,
      "STARSOLO_ALIGN " + starFull + ", STARSOLO_ALIGN_BAM_ONLY " + starBam +
      ", prior_factor " + shippedTime.STARSOLO_ALIGN?.prior_factor)
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
