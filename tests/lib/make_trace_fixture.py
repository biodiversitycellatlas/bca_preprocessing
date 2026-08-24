#!/usr/bin/env python3
"""
Generate synthetic Nextflow ``pipeline_info/`` trees for testing
``bin/resource_efficiency.py``.

The pipeline has never been run on a machine where the test suite executes, and a
real trace is several MB of cluster-specific data, so the parser and the scaling
fit are exercised against generated runs instead. Memory is drawn from a known
power law::

    peak_rss = base_memory + factor * size**exponent * lognormal(noise)

so a test can assert that the recovered exponent matches the one that was asked
for -- which is the only way to check the fit is doing arithmetic rather than
returning a plausible-looking number.

Everything is deterministic under ``--seed``.

Usage modes
-----------
Several runs of increasing size, into one outdir (what a user accumulates):
    make_trace_fixture.py --outdir /tmp/fx --runs 5 --exponent 0.85

One run per directory, for --results-root:
    make_trace_fixture.py --outdir /tmp/fx --runs 3 --separate-dirs

Example
-------
    make_trace_fixture.py --outdir /tmp/fx --runs 5 --exponent 0.85 --seed 1
    bin/resource_efficiency.py --results /tmp/fx --no-plots --json /tmp/fx.json
"""

import argparse
import datetime
import json
import math
import os
import random
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Must match the `trace.fields` list in nextflow.config: the fixture is only
# meaningful if it is shaped like what the pipeline actually writes.
TRACE_FIELDS: List[str] = [
    "task_id", "hash", "native_id", "process", "tag", "name", "status", "exit",
    "attempt", "submit", "start", "complete", "duration", "realtime", "queue",
    "cpus", "memory", "time", "%cpu", "%mem", "peak_rss", "peak_vmem", "rchar",
    "wchar", "read_bytes", "write_bytes", "error_action",
]

GB = 1024 ** 3

# A representative slice of the real pipeline. Every process name here exists in
# modules/**/main.nf (or is a genuine alias of one), so a fixture run also acts as
# a coverage check on the label mapper -- an invented name would show up as
# "unmapped" and fail the test rather than passing silently.
#
# Per entry:
#   name, fully-qualified path, requested (cpus, memory bytes, time seconds),
#   and the model (base_gb, mem_factor, time_factor_s, cpu_parallelism).
#
# Memory follows  base_gb + mem_factor * size_gb**exponent  (in GB), so the
# constant term damps the recovered exponent slightly -- exactly as a resident
# genome index does in reality. Bases are kept small relative to the scaled term
# so the embedded exponent stays recoverable within a sane tolerance.
#
# The requested values are paired with the models so that a successful task stays
# comfortably under its request: saturating the cap would flatten the curve and make
# the fit untestable. That pairing, not conf/base.config, is what fixes them.
#
# They therefore do NOT track the tier each process currently declares, and must not
# be "corrected" to -- several of these processes have been relabelled since the
# fixture was written, and matching the tiers today would put usage above request and
# destroy the property the fixture exists to test. Nothing should assert that a
# request here equals a tier value; the cases that need a real tier read it from
# conf/base.config instead (see tier_memory() in tests/checks/resource_efficiency.sh).
PROCESSES: List[Tuple[str, str, Tuple[int, int, int], Tuple[float, float, float, float]]] = [
    ("STARSOLO_ALIGN", "BCA:MAPPING:STARSOLO_ALIGN",
     (16, 128 * GB, 10 * 3600), (0.5, 1.50, 370.0, 0.72)),
    ("STARSOLO_INDEX", "BCA:MAPPING:STARSOLO_INDEX",
     (16, 128 * GB, 10 * 3600), (6.0, 0.55, 150.0, 0.55)),
    ("SALMON_INDEX", "BCA:MAPPING:SALMON_INDEX",
     (12, 64 * GB, 10 * 3600), (2.0, 0.60, 130.0, 0.60)),
    ("RM_VARBASES", "BCA:PREPROCESS:RM_VARBASES",
     (6, 36 * GB, 6 * 3600), (0.4, 0.30, 180.0, 0.80)),
    ("SAMTOOLS_VIEW_MAPPED", "BCA:POST:SAMTOOLS_VIEW_MAPPED",
     (6, 36 * GB, 6 * 3600), (0.3, 0.12, 120.0, 0.45)),
    ("SCRUBLET", "BCA:FILTER:SCRUBLET",
     (6, 36 * GB, 6 * 3600), (0.6, 0.40, 150.0, 0.30)),
    ("FASTQC", "BCA:PREPROCESS:FASTQC",
     (2, 12 * GB, 4 * 3600), (0.2, 0.09, 90.0, 0.50)),
    ("MTX_TO_10X", "BCA:POST:MTX_TO_10X",
     (2, 12 * GB, 4 * 3600), (0.2, 0.06, 60.0, 0.35)),
    # Aliased inclusion: the trace records the alias, the module declares
    # DOUBLET_FILTER. Exercises the prefix fallback in the label mapper.
    ("DOUBLET_FILTER_RAW", "BCA:FILTER:DOUBLET_FILTER_RAW",
     (2, 12 * GB, 4 * 3600), (0.2, 0.07, 70.0, 0.40)),
    # Carries no label at all, so it lands in the __default__ tier.
    ("MERGE_REF_GTF", "BCA:MAPPING:MERGE_REF_GTF",
     (1, 6 * GB, 4 * 3600), (0.1, 0.01, 25.0, 0.20)),
    ("SAVE_RUN_CONFIG", "BCA:SAVE_RUN_CONFIG",
     (1, 6 * GB, 4 * 3600), (0.05, 0.001, 5.0, 0.10)),
]


# ---------------------------------------------------------------------------
# Formatting helpers -- must produce exactly what Nextflow writes
# ---------------------------------------------------------------------------

def fmt_memory(value: float) -> str:
    """Render bytes the way Nextflow's MemoryUnit does (binary, one decimal)."""
    for unit, factor in (("TB", 1024 ** 4), ("GB", GB), ("MB", 1024 ** 2), ("KB", 1024)):
        if value >= factor:
            return f"{value / factor:.1f} {unit}"
    return f"{int(value)} B"


def fmt_duration(seconds: float) -> str:
    """Render seconds the way Nextflow's Duration does (``2h 3m 4s``)."""
    if seconds < 1:
        return f"{int(seconds * 1000)}ms"
    total = int(round(seconds))
    hours, rem = divmod(total, 3600)
    minutes, secs = divmod(rem, 60)
    parts: List[str] = []
    if hours:
        parts.append(f"{hours}h")
    if minutes:
        parts.append(f"{minutes}m")
    if secs or not parts:
        parts.append(f"{secs}s")
    return " ".join(parts)


def fmt_timestamp(moment: datetime.datetime) -> str:
    """Render a trace timestamp (``2026-08-01 09:12:33.123``)."""
    return moment.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]


# ---------------------------------------------------------------------------
# Run generation
# ---------------------------------------------------------------------------

def make_run(rng: random.Random,
             run_index: int,
             n_samples: int,
             size_gb: float,
             exponent: float,
             noise: float,
             inject_kill: bool) -> Tuple[str, List[str], List[Dict[str, str]]]:
    """Build one run's trace rows.

    Parameters
    ----------
    rng:
        Seeded random source.
    run_index:
        Position in the series, used to space the timestamps apart.
    n_samples:
        Per-sample processes emit this many tasks; reference-building ones emit one.
    size_gb:
        Nominal input size for this run. Per-sample sizes jitter around it, which
        is what puts several x values inside a single run.
    exponent:
        The power law memory follows. Tests assert this is recovered.
    noise:
        Lognormal sigma applied to memory and runtime.
    inject_kill:
        Add a killed STARSOLO_ALIGN attempt followed by a successful retry.

    Returns
    -------
    ``(run_key, trace_lines, sample_rows)``.
    """
    base_time = datetime.datetime(2026, 3, 1) + datetime.timedelta(days=7 * run_index)
    run_key = base_time.strftime("%Y-%m-%d_%H-%M-%S")

    lines: List[str] = ["\t".join(TRACE_FIELDS)]
    samples: List[Dict[str, str]] = []
    task_id = 0
    clock = base_time

    sample_sizes: List[float] = []
    for index in range(n_samples):
        # Spread samples over roughly half a decade so one run alone already has
        # an x range worth plotting.
        jitter = math.exp(rng.uniform(-0.55, 0.55))
        sample_sizes.append(size_gb * jitter)
        samples.append({
            "sample": f"run{run_index}_S{index + 1}",
            "fastq_cDNA": f"/data/run{run_index}/S{index + 1}_R2.fastq.gz",
            "fastq_BC_UMI": f"/data/run{run_index}/S{index + 1}_R1.fastq.gz",
            "expected_cells": "5000",
        })

    def emit(process: str, path: str, cpus: int, memory: int, time_limit: int,
             sample_index: Optional[int], size_bytes: float, peak_rss: float,
             realtime: float, cores: float, status: str, exit_code: int,
             attempt: int, error_action: str) -> None:
        nonlocal task_id, clock
        task_id += 1
        clock += datetime.timedelta(seconds=rng.randint(5, 120))
        start = clock
        complete = start + datetime.timedelta(seconds=realtime)
        tag = samples[sample_index]["sample"] if sample_index is not None else "-"
        row = {
            "task_id": str(task_id),
            "hash": f"{rng.randint(0, 255):02x}/{rng.randrange(16 ** 6):06x}",
            "native_id": str(4000000 + task_id + run_index * 1000),
            "process": path,
            "tag": tag,
            "name": f"{process} ({tag})" if tag != "-" else process,
            "status": status,
            "exit": str(exit_code),
            "attempt": str(attempt),
            "submit": fmt_timestamp(start - datetime.timedelta(seconds=30)),
            "start": fmt_timestamp(start),
            "complete": fmt_timestamp(complete),
            "duration": fmt_duration(realtime + 30),
            "realtime": fmt_duration(realtime),
            "queue": "genoa64",
            "cpus": str(cpus),
            "memory": fmt_memory(memory),
            "time": fmt_duration(time_limit),
            "%cpu": f"{cores * 100:.1f}%",
            "%mem": f"{100.0 * peak_rss / memory:.1f}%",
            "peak_rss": fmt_memory(peak_rss),
            "peak_vmem": fmt_memory(peak_rss * 1.4),
            "rchar": fmt_memory(size_bytes),
            "wchar": fmt_memory(size_bytes * 0.6),
            "read_bytes": fmt_memory(size_bytes * 0.9),
            "write_bytes": fmt_memory(size_bytes * 0.55),
            "error_action": error_action or "-",
        }
        lines.append("\t".join(row[field] for field in TRACE_FIELDS))
        clock = complete

    for process, path, (cpus, memory, time_limit), model in PROCESSES:
        base_gb, mem_factor, time_factor, par = model
        # Reference-building steps run once per run; everything else per sample.
        per_sample = process not in ("STARSOLO_INDEX", "SALMON_INDEX",
                                     "MERGE_REF_GTF", "SAVE_RUN_CONFIG")
        indices = range(n_samples) if per_sample else [None]

        for sample_index in indices:
            this_size = (sample_sizes[sample_index] if sample_index is not None
                         else size_gb)
            size_bytes = this_size * GB
            scaled = mem_factor * this_size ** exponent
            peak_rss = (base_gb + scaled) * GB * math.exp(rng.gauss(0, noise))
            # A COMPLETED task by definition stayed inside its request; the models
            # above are chosen so this clamp never actually binds.
            peak_rss = min(peak_rss, memory * 0.95)
            realtime = max(20.0, (30 + time_factor * this_size ** exponent)
                           * math.exp(rng.gauss(0, noise)))
            realtime = min(realtime, time_limit * 0.92)
            cores = min(cpus, max(0.1, cpus * par * math.exp(rng.gauss(0, noise / 2))))

            if (inject_kill and process == "STARSOLO_ALIGN"
                    and sample_index == 0):
                # The first attempt is killed by the cgroup, the retry succeeds with
                # task.attempt = 2 and therefore double the request.
                emit(process, path, cpus, memory, time_limit, sample_index,
                     size_bytes, memory * 0.99, realtime * 0.4, cores,
                     "FAILED", 137, 1, "RETRY")
                emit(process, path, cpus, memory * 2, time_limit * 2, sample_index,
                     size_bytes, peak_rss, realtime, cores, "COMPLETED", 0, 2, "")
                continue

            emit(process, path, cpus, memory, time_limit, sample_index,
                 size_bytes, peak_rss, realtime, cores, "COMPLETED", 0, 1, "")

    return run_key, lines, samples


def write_run(pipeline_info: str, run_key: str, lines: List[str],
              samples: List[Dict[str, str]], params: Dict[str, object]) -> None:
    """Write one run's trace, run_config, samplesheet and params JSON."""
    os.makedirs(pipeline_info, exist_ok=True)

    with open(os.path.join(pipeline_info, f"execution_trace_{run_key}.txt"),
              "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")

    with open(os.path.join(pipeline_info, f"run_config_{run_key}.txt"),
              "w", encoding="utf-8") as handle:
        for key in sorted(params):
            handle.write(f"{key} = {params[key]}\n")

    with open(os.path.join(pipeline_info, f"samplesheet_{run_key}.csv"),
              "w", encoding="utf-8", newline="") as handle:
        columns = ["sample", "fastq_cDNA", "fastq_BC_UMI", "expected_cells"]
        handle.write(",".join(columns) + "\n")
        for sample in samples:
            handle.write(",".join(sample[column] for column in columns) + "\n")

    # Deliberately offset: params_*.json carries its own timestamp, so the tool has
    # to pair it by the trace_report_suffix value inside rather than by filename.
    stamp = datetime.datetime.strptime(run_key, "%Y-%m-%d_%H-%M-%S")
    offset = (stamp + datetime.timedelta(seconds=3)).strftime("%Y-%m-%d_%H-%M-%S")
    with open(os.path.join(pipeline_info, f"params_{offset}.json"),
              "w", encoding="utf-8") as handle:
        json.dump({"trace_report_suffix": run_key,
                   "pipeline_version": "0.2.4",
                   "git_commit": "deadbeef"}, handle, indent=2)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate synthetic Nextflow pipeline_info/ trees for testing "
                    "bin/resource_efficiency.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-o", "--outdir", required=True,
                        help="Directory to create the fixture in.")
    parser.add_argument("-n", "--runs", type=int, default=5,
                        help="Number of runs to generate. Default: 5.")
    parser.add_argument("-s", "--samples", type=int, default=4,
                        help="Samples per run. Default: 4.")
    parser.add_argument("-e", "--exponent", type=float, default=0.85,
                        help="Memory scaling exponent to embed. Default: 0.85.")
    parser.add_argument("--noise", type=float, default=0.05,
                        help="Lognormal sigma on memory and runtime. Default: 0.05.")
    parser.add_argument("--min-size", type=float, default=4.0,
                        help="Input size of the smallest run, in GB. Default: 4.")
    parser.add_argument("--max-size", type=float, default=120.0,
                        help="Input size of the largest run, in GB. Default: 120.")
    parser.add_argument("--separate-dirs", action="store_true",
                        help="Give each run its own results directory, for testing "
                             "--results-root. Default: all runs share one outdir, "
                             "which is what a user accumulates in practice.")
    parser.add_argument("--no-kill", action="store_true",
                        help="Do not inject a killed task into the first run.")
    parser.add_argument("--star-limit-bamsort", default="0",
                        help="Value recorded for star_limitBAMsortRAM. Default: 0 (auto).")
    parser.add_argument("--seed", type=int, default=1,
                        help="Random seed. Default: 1.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)

    # Geometric spacing, so the runs cover a wide log range with few points --
    # which is what the fit needs and what a real series of runs looks like.
    if args.runs > 1:
        ratio = (args.max_size / args.min_size) ** (1.0 / (args.runs - 1))
        sizes = [args.min_size * ratio ** index for index in range(args.runs)]
    else:
        sizes = [args.min_size]

    for index, size in enumerate(sizes):
        run_key, lines, samples = make_run(
            rng, index, args.samples, size, args.exponent, args.noise,
            inject_kill=(index == 0 and not args.no_kill))

        results_dir = (os.path.join(args.outdir, f"outdir_run{index}")
                       if args.separate_dirs else args.outdir)
        params = {
            "protocol": "10xv3",
            "mapping_software": "starsolo",
            "outdir": results_dir,
            "trace_report_suffix": run_key,
            "star_limitBAMsortRAM": args.star_limit_bamsort,
            "star_limitGenomeGenerateRAM": "null",
        }
        write_run(os.path.join(results_dir, "pipeline_info"),
                  run_key, lines, samples, params)
        print(f"run {index}: {run_key}  ~{size:.1f} GB  "
              f"{len(lines) - 1} rows -> {results_dir}")

    print(f"\nFixture written to {args.outdir} "
          f"(exponent={args.exponent}, seed={args.seed})")


if __name__ == "__main__":
    main()
