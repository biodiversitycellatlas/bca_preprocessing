#!/usr/bin/env python3
"""
Report HPC resource efficiency across Nextflow runs and retune the scheduled
resources for the next one.

Reads the ``execution_trace_*.txt`` files that this pipeline writes into
``<outdir>/pipeline_info/`` and compares, per process, what each task actually
consumed against what the scheduler was asked for. Because the trace records
every task's own input volume (``rchar``), usage can be plotted against dataset
size -- within one run across samples of differing depth, and across runs on
differently sized datasets -- which is what makes the recommendations
extrapolate rather than merely average.

Resources in this pipeline come exclusively from the ``withLabel:`` tiers in
``conf/base.config``; no module sets ``cpus``/``memory``/``time`` directly. The
process -> label mapping is therefore recovered by scanning ``modules/**/main.nf``,
and the emitted config retunes those tiers, lifting any process that would
otherwise inflate a tier it shares into its own ``withName:`` block.

Usage modes
-----------
Report only (default; writes nothing outside --output):
    resource_efficiency.py --results /path/to/outdir

Several runs, for a scaling fit:
    resource_efficiency.py --results /path/run_A --results /path/run_B --results /path/run_C
    resource_efficiency.py --results-root /path/holding/many/outdirs

Emit a tuned config and wire it into the next run:
    resource_efficiency.py --results /path/to/outdir --apply
    # writes conf/resources_tuned.config, which submit_nextflow.sh picks up

Extrapolate to a dataset larger than anything observed:
    resource_efficiency.py --results-root /path/runs --target-size 400

Text only, no matplotlib needed:
    resource_efficiency.py --results /path/to/outdir --no-plots

Notes
-----
The analysis core is standard library only, so it runs on a bare login-node
``python3``. matplotlib is needed only for the figures; without it, pass
``--no-plots`` and the table, TSV and config emission still work.
"""

import argparse
import base64
import csv
import datetime
import glob
import json
import math
import os
import re
import sys
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Tiers in conf/base.config that carry cpus/memory/time. Labels outside this set
# (error_ignore, error_optional, error_retry, process_gpu) change error handling
# or accelerators, not resource size, so they are recorded but never retuned.
#
# This set is DISCOVERED from conf/base.config by ``adopt_resource_labels()``, not
# maintained by hand: the hardcoded list below is only the fallback for when
# base.config cannot be read. A hand-maintained list silently rots -- when the *2
# tiers were added, every process carrying one was mapped to DEFAULT_TIER, and the
# resulting 6 GB denominator corrupted every efficiency figure and every emitted
# withLabel block until it was noticed.
FALLBACK_RESOURCE_LABELS: List[str] = [
    "process_single",
    "process_low",
    "process_medium",
    "process_high",
    "process_high_memory",
    "process_long",
]

RESOURCE_LABELS: List[str] = list(FALLBACK_RESOURCE_LABELS)

# Every withLabel block conf/base.config declares, resource-bearing or not. Used to
# tell "declared a label this tool does not know" apart from "declared no label";
# empty until base.config has been read.
KNOWN_LABELS: Set[str] = set()

# Pseudo-tier for processes carrying no label at all (MERGE_REF_FASTA,
# MERGE_REF_GTF). They inherit the bare `process {}` block, which is shared with
# every unlabelled process, so they are never retuned as a tier -- only ever as
# their own withName block.
DEFAULT_TIER = "__default__"

# Trace statuses. CACHED rows come from `-resume` and carry '-' for every usage
# field; counting them would both double-count the task and drag the stats to
# zero. FAILED is kept aside rather than dropped -- it is the strongest
# under-provisioning signal in the file.
STATUS_COMPLETED = "COMPLETED"
STATUS_FAILED = "FAILED"

# Exit codes that mean "the scheduler killed this": 137 = SIGKILL (OOM killer or
# cgroup limit), 140 = SLURM walltime on some sites, 143 = SIGTERM (graceful
# walltime), 125/104 appear on cgroup-v2 memory kills.
KILLED_EXIT_CODES = {104, 125, 137, 140, 143}

# Safety factors applied to the observed statistic before recommending.
DEFAULT_SAFETY_MEMORY = 1.25   # /proc polling can miss a short allocation spike
DEFAULT_SAFETY_TIME = 1.50     # walltime overrun kills the task; slack is cheap

# Floors and quantisation of the emitted numbers. Memory is quantised coarsely
# only once it is large: rounding a 200 MB process up to the next 4 GB would waste
# more than it uses, so small requests step in 1 GB instead.
MIN_MEMORY_BYTES = 1 * 1024 ** 3
MEMORY_STEP_BYTES = 4 * 1024 ** 3
MEMORY_FINE_STEP_BYTES = 1 * 1024 ** 3
MEMORY_FINE_BELOW_BYTES = 16 * 1024 ** 3
# Walltime steps in half hours while it is short and in whole hours once it is
# long, so a long request reads as "13.h" rather than "750.m".
MIN_TIME_SECONDS = 30 * 60
TIME_STEP_SECONDS = 30 * 60
TIME_COARSE_STEP_SECONDS = 3600
TIME_COARSE_ABOVE_SECONDS = 4 * 3600

# A log-log fit is only believed when it clears all four of these.
MIN_FIT_POINTS = 4      # distinct tasks
MIN_FIT_RUNS = 2        # so the slope is not one run's sample spread
MIN_FIT_R2 = 0.7
MIN_FIT_LOG_SPAN = 0.3  # x range spans at least ~2x

# A tier member needing more than this multiple of the tier median is lifted out
# into its own withName block instead of inflating the shared tier.
OUTLIER_FACTOR = 1.5
# Below this many members the median is meaningless, so nothing is called an outlier.
MIN_MEMBERS_FOR_OUTLIER = 3

# Tier values within this fraction of the current declaration are left alone.
NO_CHANGE_TOLERANCE = 0.10

# Processes whose tool imposes its own memory ceiling, independent of what the
# scheduler allocated. Raising the Nextflow `memory` directive for these does
# nothing on its own while the cap stays pinned below it, and lowering the
# directive below a pinned cap gets the task killed by the cgroup -- so a
# recommendation here is only actionable once the cap is checked too.
#   {process: (param name, human note)}
TOOL_MEMORY_CAPS: Dict[str, Tuple[str, str]] = {
    "STARSOLO_ALIGN": (
        "star_limitBAMsortRAM",
        "STAR's BAM sort buffer. Left at 0 STAR sets it to the genome index size, "
        "which is unrelated to the allocation -- the classic 'not enough memory for "
        "BAM sorting' failure on an otherwise idle 128 GB job.",
    ),
    "STARSOLO_INDEX": (
        "star_limitGenomeGenerateRAM",
        "STAR's ceiling for suffix-array construction. A pinned value neither grows "
        "with a larger tier nor shrinks with a smaller one.",
    ),
}

# Values of a cap param that mean "derive it from task.memory" rather than a pin.
CAP_AUTO_TOKENS = {"", "0", "null", "none", "auto"}

BYTE_UNITS: Dict[str, int] = {
    "B": 1,
    "KB": 1024,
    "MB": 1024 ** 2,
    "GB": 1024 ** 3,
    "TB": 1024 ** 4,
    "PB": 1024 ** 5,
    # Nextflow's MemoryUnit prints KB/MB/GB but means KiB/MiB/GiB; accept the
    # explicit spellings too in case a trace was produced by other tooling.
    "KIB": 1024,
    "MIB": 1024 ** 2,
    "GIB": 1024 ** 3,
    "TIB": 1024 ** 4,
}

TIME_UNITS: Dict[str, float] = {
    "ms": 0.001,
    "s": 1.0,
    "m": 60.0,
    "h": 3600.0,
    "d": 86400.0,
}

# Values Nextflow writes when a field was never measured.
NULL_TOKENS = {"", "-", "na", "n/a", "null", "none"}

_DURATION_TOKEN_RE = re.compile(r"([0-9]*\.?[0-9]+)\s*(ms|s|m|h|d)", re.IGNORECASE)
_PROCESS_DECL_RE = re.compile(r"^\s*process\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{")
_LABEL_RE = re.compile(r"""^\s*label\s*:?\s*['"]([A-Za-z0-9_]+)['"]""")
_SECTION_RE = re.compile(r"^\s*(input|output|script|shell|exec|when|stub)\s*:")


# ---------------------------------------------------------------------------
# Unit parsing
#
# Every parser returns None rather than raising, so one malformed cell degrades
# a single metric instead of aborting the report.
# ---------------------------------------------------------------------------

def _is_null(text: Optional[str]) -> bool:
    """True when *text* is one of Nextflow's not-measured placeholders."""
    return text is None or text.strip().lower() in NULL_TOKENS


def parse_bytes(text: Optional[str]) -> Optional[float]:
    """Parse a Nextflow memory cell into bytes.

    Handles ``"1.2 GB"``, ``"512 MB"``, ``"6GB"`` and, in ``trace.raw`` mode, a
    bare integer count of bytes. Multipliers are **base 1024**: Nextflow's
    ``MemoryUnit`` prints ``GB`` but means ``GiB``, and treating it as 10^9
    understates every figure by about 7%.

    Parameters
    ----------
    text:
        Raw cell contents.

    Returns
    -------
    Bytes as a float, or ``None`` when the cell holds no measurement.
    """
    if _is_null(text):
        return None
    raw = text.strip()
    match = re.fullmatch(r"([0-9]*\.?[0-9]+)\s*([A-Za-z]*)", raw)
    if not match:
        return None
    value, unit = float(match.group(1)), match.group(2).upper()
    if not unit:
        return value
    return value * BYTE_UNITS.get(unit, 1) if unit in BYTE_UNITS else None


def parse_duration(text: Optional[str], raw_mode: bool = False) -> Optional[float]:
    """Parse a Nextflow duration cell into seconds.

    Handles the compound form (``"2h 3m 4s"``, ``"1d 4h"``), a single unit
    (``"45.2s"``, ``"250ms"``) and a bare number. A bare number is **milliseconds**
    -- that is what ``trace.raw = true`` emits, and Nextflow never writes a
    unitless seconds value.

    Parameters
    ----------
    text:
        Raw cell contents.
    raw_mode:
        Unused for the decision (a bare number is always ms) but kept explicit so
        callers document which trace flavour they believe they are reading.

    Returns
    -------
    Seconds as a float, or ``None``.
    """
    if _is_null(text):
        return None
    raw = text.strip()
    tokens = _DURATION_TOKEN_RE.findall(raw)
    if tokens:
        # Reject a partial match such as "12x3h": every character must be consumed
        # by the tokens plus separators.
        if re.fullmatch(r"(\s*[0-9]*\.?[0-9]+\s*(?:ms|s|m|h|d)\s*)+", raw, re.IGNORECASE):
            return sum(float(v) * TIME_UNITS[u.lower()] for v, u in tokens)
        return None
    try:
        return float(raw) / 1000.0
    except ValueError:
        return None


def parse_percent(text: Optional[str]) -> Optional[float]:
    """Parse a ``%cpu``/``%mem`` cell.

    The value may exceed 100: ``%cpu`` is summed across cores, so a task using
    1.4 of its 6 allocated cores reads ``140.0%``.
    """
    if _is_null(text):
        return None
    try:
        return float(text.strip().rstrip("%"))
    except ValueError:
        return None


def parse_int(text: Optional[str]) -> Optional[int]:
    """Parse an integer cell, tolerating floats and placeholders."""
    if _is_null(text):
        return None
    try:
        return int(float(text.strip()))
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def fmt_bytes(value: Optional[float]) -> str:
    """Render bytes as a compact binary-unit string."""
    if value is None:
        return "n/a"
    for unit, factor in (("TB", 1024 ** 4), ("GB", 1024 ** 3), ("MB", 1024 ** 2), ("KB", 1024)):
        if value >= factor:
            return f"{value / factor:.1f} {unit}"
    return f"{value:.0f} B"


def fmt_duration(value: Optional[float]) -> str:
    """Render seconds as ``1h 5m`` / ``45s``."""
    if value is None:
        return "n/a"
    seconds = int(round(value))
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes}m"
    if minutes:
        return f"{minutes}m {secs}s"
    return f"{secs}s"


def fmt_pct(value: Optional[float]) -> str:
    """Render a 0-1 ratio as a percentage."""
    return "n/a" if value is None else f"{value * 100:.0f}%"


def groovy_memory(value: float) -> str:
    """Render bytes as a Groovy memory literal (``28.GB``)."""
    gb = value / 1024 ** 3
    if gb >= 1 and abs(gb - round(gb)) < 1e-9:
        return f"{int(round(gb))}.GB"
    mb = value / 1024 ** 2
    return f"{int(round(mb))}.MB"


def groovy_time(value: float) -> str:
    """Render seconds as a Groovy duration literal (``6.h`` / ``90.m``)."""
    if value % 3600 == 0:
        return f"{int(value // 3600)}.h"
    return f"{int(round(value / 60))}.m"


# ---------------------------------------------------------------------------
# Small statistics helpers
#
# Deliberately dependency-free: the analysis has to run on whatever python3 a
# login node happens to offer, and none of this is heavy enough to justify numpy.
# ---------------------------------------------------------------------------

def percentile(values: Sequence[float], q: float) -> Optional[float]:
    """Linear-interpolated percentile of *values* at *q* (0-100)."""
    data = sorted(v for v in values if v is not None)
    if not data:
        return None
    if len(data) == 1:
        return data[0]
    pos = (len(data) - 1) * q / 100.0
    lower, upper = math.floor(pos), math.ceil(pos)
    if lower == upper:
        return data[int(pos)]
    return data[lower] * (upper - pos) + data[upper] * (pos - lower)


def median(values: Sequence[float]) -> Optional[float]:
    """Median of *values*, or ``None`` when empty."""
    return percentile(values, 50)


class FitResult:
    """Outcome of a log-log ordinary-least-squares fit ``y = a * x**b``."""

    def __init__(self, slope: float, intercept: float, r2: float,
                 resid_std: float, n_points: int, n_runs: int, log_span: float) -> None:
        self.slope = slope
        self.intercept = intercept
        self.r2 = r2
        self.resid_std = resid_std
        self.n_points = n_points
        self.n_runs = n_runs
        self.log_span = log_span

    @property
    def trusted(self) -> bool:
        """Whether the fit clears every quality gate.

        All four matter: too few points overfits, a single run measures sample
        spread rather than scaling, a poor r2 means the relationship is not a
        power law, and a narrow x range makes the slope an extrapolation from
        noise.
        """
        return (self.n_points >= MIN_FIT_POINTS
                and self.n_runs >= MIN_FIT_RUNS
                and self.r2 >= MIN_FIT_R2
                and self.log_span >= MIN_FIT_LOG_SPAN)

    def predict(self, x_value: float, sigma: float = 2.0) -> float:
        """Predict y at *x_value*, padded by *sigma* residual standard deviations."""
        log_pred = self.intercept + self.slope * math.log10(x_value)
        return 10 ** (log_pred + sigma * self.resid_std)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "exponent": self.slope,
            "intercept": self.intercept,
            "r2": self.r2,
            "resid_std": self.resid_std,
            "n_points": self.n_points,
            "n_runs": self.n_runs,
            "log_span": self.log_span,
            "trusted": self.trusted,
        }


def fit_loglog(points: Sequence[Tuple[float, float, str]]) -> Optional[FitResult]:
    """Fit ``y = a * x**b`` by OLS on log10 axes.

    Parameters
    ----------
    points:
        ``(x, y, run_key)`` triples. Non-positive x or y are dropped, since the
        fit is on log axes.

    Returns
    -------
    A :class:`FitResult`, or ``None`` when there is nothing to fit.
    """
    usable = [(math.log10(x), math.log10(y), key)
              for x, y, key in points if x and y and x > 0 and y > 0]
    if len(usable) < 2:
        return None

    xs = [p[0] for p in usable]
    ys = [p[1] for p in usable]
    n = len(usable)
    mean_x, mean_y = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    syy = sum((y - mean_y) ** 2 for y in ys)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    if sxx <= 0:
        return None

    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    r2 = (sxy * sxy) / (sxx * syy) if syy > 0 else 1.0
    residuals = [y - (slope * x + intercept) for x, y in zip(xs, ys)]
    resid_std = math.sqrt(sum(r * r for r in residuals) / max(1, n - 2))

    return FitResult(
        slope=slope,
        intercept=intercept,
        r2=r2,
        resid_std=resid_std,
        n_points=n,
        n_runs=len({p[2] for p in usable}),
        log_span=max(xs) - min(xs),
    )


# ---------------------------------------------------------------------------
# Label mapping
#
# The trace records no label, so the process -> tier mapping has to be rebuilt
# from the module sources and conf/base.config.
# ---------------------------------------------------------------------------

class ProcessLabels:
    """Labels declared by one ``process`` block in a module file."""

    def __init__(self, name: str, labels: List[str], source: str) -> None:
        self.name = name
        self.labels = labels
        self.source = source

    @property
    def tier(self) -> str:
        """The resource-bearing label, or ``DEFAULT_TIER`` when there is none."""
        for label in self.labels:
            if label in RESOURCE_LABELS:
                return label
        return DEFAULT_TIER

    @property
    def modifiers(self) -> List[str]:
        """Labels that are not resource tiers (error_*, process_gpu)."""
        return [label for label in self.labels if label not in RESOURCE_LABELS]

    @property
    def unrecognised(self) -> List[str]:
        """Labels this process declares that conf/base.config does not define at all.

        ``tier`` cannot express this on its own: it answers ``DEFAULT_TIER`` both for
        a process that declares no label and for one that declares a label this tool
        has never heard of. Those are very different -- the second means the report
        is measuring against the wrong tier -- so they are separated here, and
        ``tests/checks/resource_efficiency.sh::label_coverage`` asserts this is empty.

        Empty while ``KNOWN_LABELS`` is unpopulated, so that callers which never read
        conf/base.config do not report every label as unknown.
        """
        if not KNOWN_LABELS:
            return []
        return [label for label in self.labels if label not in KNOWN_LABELS]


def scan_module_labels(pipeline_dir: str) -> Dict[str, ProcessLabels]:
    """Map every process name declared under ``modules/`` to its labels.

    A module may declare several labels (``CELLBENDER`` is ``process_high`` plus
    ``process_gpu``; ``SCRUBLET`` is ``process_medium`` plus ``error_optional``),
    so all are collected and :attr:`ProcessLabels.tier` picks the resource one.

    Parameters
    ----------
    pipeline_dir:
        Pipeline checkout containing ``modules/``.

    Returns
    -------
    ``{process_name: ProcessLabels}``.
    """
    found: Dict[str, ProcessLabels] = {}
    pattern = os.path.join(pipeline_dir, "modules", "**", "*.nf")

    for path in glob.glob(pattern, recursive=True):
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                lines = handle.readlines()
        except OSError as exc:
            sys.stderr.write(f"Warning: could not read {path}: {exc}\n")
            continue

        current: Optional[str] = None
        labels: List[str] = []
        for line in lines:
            decl = _PROCESS_DECL_RE.match(line)
            if decl:
                if current:
                    found[current] = ProcessLabels(current, labels, path)
                current, labels = decl.group(1), []
                continue
            if current is None:
                continue
            # Labels are declared in the directive block, before input:/script:.
            if _SECTION_RE.match(line):
                found[current] = ProcessLabels(current, labels, path)
                current, labels = None, []
                continue
            label = _LABEL_RE.match(line)
            if label:
                labels.append(label.group(1))
        if current:
            found[current] = ProcessLabels(current, labels, path)

    return found


class TierSpec:
    """Resources a ``withLabel:`` block declares in conf/base.config."""

    def __init__(self, label: str) -> None:
        self.label = label
        self.cpus: Optional[float] = None
        self.memory: Optional[float] = None
        self.time: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return {"label": self.label, "cpus": self.cpus,
                "memory": self.memory, "time": self.time}


def _groovy_value(expression: str) -> str:
    """Reduce ``{ 36.GB * task.attempt }`` to the first-attempt literal ``36.GB``.

    The closing brace has to be stripped explicitly: a directive without a
    multiplier (``cpus = { 1 }``, as process_single declares) leaves nothing for
    the ``*`` split to remove.
    """
    return expression.strip().lstrip("{").split("*")[0].strip().rstrip("}").strip()


def _parse_groovy_memory(expression: str) -> Optional[float]:
    """Parse a Groovy memory literal (``36.GB``, ``512.MB``) into bytes."""
    value = _groovy_value(expression)
    match = re.match(r"^([0-9]*\.?[0-9]+)\s*\.?\s*([A-Za-z]+)$", value)
    if match:
        return parse_bytes(f"{match.group(1)} {match.group(2)}")
    return parse_bytes(value)


def _parse_groovy_time(expression: str) -> Optional[float]:
    """Parse a Groovy duration literal (``6.h``, ``20.h``, ``90.m``) into seconds."""
    value = _groovy_value(expression)
    match = re.match(r"^([0-9]*\.?[0-9]+)\s*\.\s*([a-zA-Z]+)$", value)
    if match:
        return parse_duration(f"{match.group(1)}{match.group(2)}")
    return parse_duration(value)


def _extract_block(text: str, open_index: int) -> str:
    """Return the contents of the brace-delimited block starting at *open_index*.

    A regex cannot do this: every directive in a tier is itself a closure, so
    ``[^}]*`` stops at the ``}`` of ``{ 6 * task.attempt }`` and the memory and
    time lines are never seen. Depth counting is the only correct reading.
    """
    depth = 0
    for position in range(open_index, len(text)):
        if text[position] == "{":
            depth += 1
        elif text[position] == "}":
            depth -= 1
            if depth == 0:
                return text[open_index + 1:position]
    return text[open_index + 1:]


def parse_base_config(path: str) -> Dict[str, TierSpec]:
    """Extract each tier's declared cpus/memory/time from conf/base.config.

    Values look like ``{ 36.GB * task.attempt }``; the ``task.attempt`` multiplier
    is stripped, leaving the first-attempt request. The unlabelled ``process {}``
    defaults are recorded under :data:`DEFAULT_TIER`.

    Parameters
    ----------
    path:
        Path to ``conf/base.config``.

    Returns
    -------
    ``{label: TierSpec}``. Empty when the file is unreadable.
    """
    tiers: Dict[str, TierSpec] = {}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    except OSError as exc:
        sys.stderr.write(f"Warning: could not read {path}: {exc}\n")
        return tiers

    def read_directives(block: str, spec: TierSpec) -> None:
        # Only top-level directives of this block: a nested withLabel inside the
        # process{} block must not leak into the defaults.
        for key in ("cpus", "memory", "time"):
            match = re.search(r"^\s*" + key + r"\s*=\s*(.+)$", block, re.MULTILINE)
            if not match:
                continue
            expression = match.group(1)
            if key == "cpus":
                try:
                    spec.cpus = float(_groovy_value(expression))
                except ValueError:
                    pass
            elif key == "memory":
                spec.memory = _parse_groovy_memory(expression)
            else:
                spec.time = _parse_groovy_time(expression)

    # The bare `process { ... }` defaults are the directives before the first
    # withLabel block; everything after belongs to a tier.
    process_match = re.search(r"^\s*process\s*\{", text, re.MULTILINE)
    if process_match:
        body = _extract_block(text, process_match.end() - 1)
        head = body.split("withLabel")[0]
        default_spec = TierSpec(DEFAULT_TIER)
        read_directives(head, default_spec)
        tiers[DEFAULT_TIER] = default_spec

    for match in re.finditer(r"withLabel\s*:\s*([A-Za-z0-9_]+)\s*\{", text):
        label = match.group(1)
        spec = TierSpec(label)
        read_directives(_extract_block(text, match.end() - 1), spec)
        tiers[label] = spec

    return tiers


def adopt_resource_labels(tiers: Dict[str, TierSpec]) -> List[str]:
    """Rebind :data:`RESOURCE_LABELS` and :data:`KNOWN_LABELS` from parsed *tiers*.

    A tier is resource-bearing exactly when it declares at least one of
    cpus/memory/time; the behaviour-only labels (``error_ignore``,
    ``error_optional``, ``error_retry``, ``process_gpu``) declare none of the three
    and are therefore excluded, which is the same distinction the hand-maintained
    list used to encode -- only now it cannot fall out of date.

    Declaration order in ``conf/base.config`` is preserved, since that is the order
    the tier sections of the report are rendered in.

    Falls back to :data:`FALLBACK_RESOURCE_LABELS` when *tiers* is empty (an
    unreadable base.config), so the tool still runs against a trace alone.
    """
    global RESOURCE_LABELS, KNOWN_LABELS

    if not tiers:
        RESOURCE_LABELS = list(FALLBACK_RESOURCE_LABELS)
        KNOWN_LABELS = set(FALLBACK_RESOURCE_LABELS)
        return RESOURCE_LABELS

    RESOURCE_LABELS = [
        name for name, spec in tiers.items()
        if name != DEFAULT_TIER
        and (spec.cpus is not None or spec.memory is not None or spec.time is not None)
    ]
    KNOWN_LABELS = {name for name in tiers if name != DEFAULT_TIER}
    return RESOURCE_LABELS


def simple_process_name(process_field: str, name_field: str) -> str:
    """Reduce a trace row's process columns to a bare process name.

    The ``process`` column carries the fully qualified path
    (``NFCORE_BCA:BCA:MAPPING:STARSOLO_ALIGN``); ``name`` carries the same plus
    the tag (``STARSOLO_ALIGN (sample_A)``). Either is acceptable, so the more
    reliable one is preferred and the other used as a fallback.
    """
    candidate = process_field if not _is_null(process_field) else name_field
    if _is_null(candidate):
        return "UNKNOWN"
    candidate = re.sub(r"\s*\(.*\)\s*$", "", candidate.strip())
    return candidate.split(":")[-1].strip()


def resolve_label(process_name: str,
                  module_labels: Dict[str, ProcessLabels]) -> Tuple[Optional[str], Optional[str]]:
    """Map a trace process name to its declared module and tier.

    Exact match first. Failing that, the longest known process name that is a
    prefix of *process_name*, which recovers the aliased inclusions -- the trace
    records ``DOUBLET_FILTER_RAW`` and ``MERGE_REF_GTF_GENEEXT``, but the module
    files only ever declare ``DOUBLET_FILTER`` and ``MERGE_REF_GTF``.

    Returns
    -------
    ``(module_process_name, tier)``, or ``(None, None)`` when unresolvable, in
    which case the caller must surface it rather than silently drop the task.
    """
    if process_name in module_labels:
        return process_name, module_labels[process_name].tier

    candidates = [known for known in module_labels if process_name.startswith(known)]
    if candidates:
        best = max(candidates, key=len)
        return best, module_labels[best].tier

    return None, None


# ---------------------------------------------------------------------------
# Run discovery
# ---------------------------------------------------------------------------

class RunContext:
    """One previous pipeline run, identified by its trace file's timestamp suffix."""

    def __init__(self, key: str, trace_path: str, results_dir: str) -> None:
        self.key = key
        self.trace_path = trace_path
        self.results_dir = results_dir
        self.run_config_path: Optional[str] = None
        self.samplesheet_path: Optional[str] = None
        self.params: Dict[str, Any] = {}
        self.samples: List[Dict[str, str]] = []

    @property
    def label(self) -> str:
        """Short human-facing name for legends and tables."""
        return self.key

    def describe(self) -> str:
        protocol = self.params.get("protocol") or "?"
        mapper = self.params.get("mapping_software") or "?"
        return f"{self.key}  protocol={protocol}  mapper={mapper}  samples={len(self.samples)}"


def _read_run_config(path: str) -> Dict[str, str]:
    """Parse a ``run_config_<suffix>.txt`` of sorted ``key = value`` lines."""
    params: Dict[str, str] = {}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if "=" not in line:
                    continue
                key, _, value = line.partition("=")
                params[key.strip()] = value.strip()
    except OSError as exc:
        sys.stderr.write(f"Warning: could not read {path}: {exc}\n")
    return params


def _read_samplesheet(path: str) -> List[Dict[str, str]]:
    """Parse a ``samplesheet_<suffix>.csv`` into row dicts."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
            return [row for row in csv.DictReader(handle)]
    except OSError as exc:
        sys.stderr.write(f"Warning: could not read {path}: {exc}\n")
        return []


def _pair_by_suffix(pipeline_info: str, prefix: str, extension: str,
                    run_key: str) -> Optional[str]:
    """Find the ``<prefix><suffix><extension>`` file belonging to *run_key*.

    An exact suffix match is preferred, but it cannot be relied on. ``-resume``
    re-evaluates ``params.trace_report_suffix``, so the resumed run writes a new
    ``execution_trace_<newkey>.txt`` while ``SAVE_RUN_CONFIG`` is cached and never
    republishes -- leaving the run config under the *original* timestamp. Falling
    back to the nearest earlier suffix recovers the file that was actually current
    when the trace was written; without it a resumed run reports no samples and no
    protocol at all.

    Parameters
    ----------
    pipeline_info:
        The run's ``pipeline_info/`` directory.
    prefix, extension:
        Filename parts around the timestamp, e.g. ``run_config_`` and ``.txt``.
    run_key:
        The trace's timestamp suffix.

    Returns
    -------
    Best-matching path, or ``None`` when no candidate exists.
    """
    exact = os.path.join(pipeline_info, f"{prefix}{run_key}{extension}")
    if os.path.exists(exact):
        return exact

    candidates: List[Tuple[str, str]] = []
    for path in glob.glob(os.path.join(pipeline_info, f"{prefix}*{extension}")):
        name = os.path.basename(path)
        suffix = name[len(prefix):-len(extension)] if extension else name[len(prefix):]
        candidates.append((suffix, path))
    if not candidates:
        return None

    # The suffix is 'YYYY-MM-DD_HH-MM-SS', so lexicographic order is chronological.
    candidates.sort()
    earlier = [path for suffix, path in candidates if suffix <= run_key]
    return earlier[-1] if earlier else candidates[0][1]


def _pair_params_json(pipeline_info: str, run_key: str) -> Dict[str, Any]:
    """Find the ``params_*.json`` belonging to *run_key*.

    These files are timestamped by their own ``new Date()`` call, not by
    ``params.trace_report_suffix``, so their filename is off by seconds from the
    trace's. They do however *contain* ``trace_report_suffix``, so the pairing is
    done on that value rather than on filename proximity.
    """
    for path in sorted(glob.glob(os.path.join(pipeline_info, "params_*.json"))):
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                payload = json.load(handle)
        except (OSError, ValueError):
            continue
        for container in (payload, payload.get("params", {}) if isinstance(payload, dict) else {}):
            if isinstance(container, dict) and container.get("trace_report_suffix") == run_key:
                return payload
    return {}


def discover_runs(results_dirs: Sequence[str], results_root: Optional[str]) -> List[RunContext]:
    """Collect every run whose trace is reachable from the given directories.

    Parameters
    ----------
    results_dirs:
        Explicit ``params.outdir`` paths.
    results_root:
        A directory whose immediate children are treated as result directories.

    Returns
    -------
    Runs sorted by key (i.e. chronologically, the suffix being a timestamp).
    """
    candidates: List[str] = list(results_dirs)
    if results_root:
        for entry in sorted(os.listdir(results_root)):
            full = os.path.join(results_root, entry)
            if os.path.isdir(full):
                candidates.append(full)

    runs: Dict[str, RunContext] = {}
    for results_dir in candidates:
        pipeline_info = os.path.join(results_dir, "pipeline_info")
        if not os.path.isdir(pipeline_info):
            # Tolerate being pointed straight at a pipeline_info directory.
            if os.path.basename(os.path.normpath(results_dir)) == "pipeline_info":
                pipeline_info = results_dir
                results_dir = os.path.dirname(os.path.normpath(results_dir))
            else:
                sys.stderr.write(f"Warning: no pipeline_info/ under {results_dir}, skipping\n")
                continue

        traces = sorted(glob.glob(os.path.join(pipeline_info, "execution_trace_*.txt")))
        if not traces:
            sys.stderr.write(
                f"Warning: no execution_trace_*.txt in {pipeline_info}.\n"
                f"         Runs from before the trace block was restored in nextflow.config\n"
                f"         produce no trace; re-run the pipeline to collect one.\n")
            continue

        for trace_path in traces:
            base = os.path.basename(trace_path)
            key = base[len("execution_trace_"):-len(".txt")]
            if key in runs:
                continue
            run = RunContext(key, trace_path, results_dir)

            run_config = _pair_by_suffix(pipeline_info, "run_config_", ".txt", key)
            if run_config:
                run.run_config_path = run_config
                run.params = dict(_read_run_config(run_config))

            samplesheet = _pair_by_suffix(pipeline_info, "samplesheet_", ".csv", key)
            if samplesheet:
                run.samplesheet_path = samplesheet
                run.samples = _read_samplesheet(samplesheet)

            payload = _pair_params_json(pipeline_info, key)
            for field in ("pipeline_version", "version", "git_commit", "commit"):
                if isinstance(payload, dict) and field in payload:
                    run.params.setdefault(field, str(payload[field]))

            runs[key] = run

    return [runs[key] for key in sorted(runs)]


# ---------------------------------------------------------------------------
# Trace parsing
# ---------------------------------------------------------------------------

class TaskRecord:
    """One row of an execution trace, normalised to bytes and seconds."""

    __slots__ = (
        "run_key", "task_id", "process", "name", "tag", "status", "exit_code",
        "attempt", "realtime", "duration", "req_cpus", "req_memory", "req_time",
        "pct_cpu", "pct_mem", "peak_rss", "peak_vmem", "rchar", "wchar",
        "native_id", "queue", "error_action", "module_process", "tier",
    )

    def __init__(self) -> None:
        self.run_key = ""
        self.task_id = ""
        self.process = ""
        self.name = ""
        self.tag = ""
        self.status = ""
        self.exit_code: Optional[int] = None
        self.attempt = 1
        self.realtime: Optional[float] = None
        self.duration: Optional[float] = None
        self.req_cpus: Optional[float] = None
        self.req_memory: Optional[float] = None
        self.req_time: Optional[float] = None
        self.pct_cpu: Optional[float] = None
        self.pct_mem: Optional[float] = None
        self.peak_rss: Optional[float] = None
        self.peak_vmem: Optional[float] = None
        self.rchar: Optional[float] = None
        self.wchar: Optional[float] = None
        self.native_id = ""
        self.queue = ""
        self.error_action = ""
        self.module_process: Optional[str] = None
        self.tier: Optional[str] = None

    @property
    def killed(self) -> bool:
        """Whether this task looks scheduler-killed (OOM or walltime)."""
        return (self.exit_code in KILLED_EXIT_CODES
                or self.error_action.upper() == "RETRY")

    @property
    def cpu_cores_used(self) -> Optional[float]:
        """Cores actually kept busy, derived from the summed-across-cores %cpu."""
        return None if self.pct_cpu is None else self.pct_cpu / 100.0

    def memory_efficiency(self) -> Optional[float]:
        if self.peak_rss is None or not self.req_memory:
            return None
        return self.peak_rss / self.req_memory

    def cpu_efficiency(self) -> Optional[float]:
        if self.pct_cpu is None or not self.req_cpus:
            return None
        return (self.pct_cpu / 100.0) / self.req_cpus

    def time_efficiency(self) -> Optional[float]:
        if self.realtime is None or not self.req_time:
            return None
        return self.realtime / self.req_time


class TraceParseResult:
    """Everything one trace file yielded, including what could not be read."""

    def __init__(self, run_key: str) -> None:
        self.run_key = run_key
        self.completed: List[TaskRecord] = []
        self.failures: List[TaskRecord] = []
        self.n_cached = 0
        self.n_other = 0
        self.missing_fields: List[str] = []
        self.raw_mode = False


def parse_trace(path: str, run_key: str) -> TraceParseResult:
    """Read one ``execution_trace_*.txt`` into normalised task records.

    Columns are located by name from the header, never by position, because a
    trace produced with a different ``trace.fields`` list (or by a collaborator's
    pipeline) reorders them freely. Fields the file does not carry are reported in
    :attr:`TraceParseResult.missing_fields` and degrade the affected metric only.

    Parameters
    ----------
    path:
        Trace file to read.
    run_key:
        Timestamp suffix identifying the run, stamped onto every record.

    Returns
    -------
    A :class:`TraceParseResult`.
    """
    result = TraceParseResult(run_key)
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    except OSError as exc:
        sys.stderr.write(f"Warning: could not read {path}: {exc}\n")
        return result

    if not rows:
        return result

    header = [cell.strip() for cell in rows[0]]
    index = {field: position for position, field in enumerate(header)}

    for expected in ("cpus", "memory", "time", "attempt", "peak_rss", "realtime", "rchar"):
        if expected not in index:
            result.missing_fields.append(expected)

    def cell(row: List[str], field: str) -> Optional[str]:
        position = index.get(field)
        if position is None or position >= len(row):
            return None
        return row[position]

    # trace.raw = true writes every memory as a bare byte count. Detecting it up
    # front means durations in the same file are read as milliseconds.
    memory_cells = [cell(row, "memory") for row in rows[1:]]
    populated = [value for value in memory_cells if not _is_null(value)]
    if populated and all(re.fullmatch(r"\d+", value.strip()) for value in populated):
        result.raw_mode = True

    for row in rows[1:]:
        status = (cell(row, "status") or "").strip().upper()

        if status == "CACHED":
            # A -resume writes these with '-' in every usage column; counting them
            # would double-count the task and drag every statistic toward zero.
            result.n_cached += 1
            continue
        if status not in (STATUS_COMPLETED, STATUS_FAILED):
            result.n_other += 1
            continue

        record = TaskRecord()
        record.run_key = run_key
        record.task_id = (cell(row, "task_id") or "").strip()
        record.process = (cell(row, "process") or "").strip()
        record.name = (cell(row, "name") or "").strip()
        raw_tag = cell(row, "tag")
        record.tag = "" if _is_null(raw_tag) else raw_tag.strip()
        record.status = status
        record.exit_code = parse_int(cell(row, "exit"))
        record.attempt = parse_int(cell(row, "attempt")) or 1
        record.realtime = parse_duration(cell(row, "realtime"), result.raw_mode)
        record.duration = parse_duration(cell(row, "duration"), result.raw_mode)
        record.req_cpus = parse_int(cell(row, "cpus"))
        record.req_memory = parse_bytes(cell(row, "memory"))
        record.req_time = parse_duration(cell(row, "time"), result.raw_mode)
        record.pct_cpu = parse_percent(cell(row, "%cpu"))
        record.pct_mem = parse_percent(cell(row, "%mem"))
        record.peak_rss = parse_bytes(cell(row, "peak_rss"))
        record.peak_vmem = parse_bytes(cell(row, "peak_vmem"))
        record.rchar = parse_bytes(cell(row, "rchar"))
        record.wchar = parse_bytes(cell(row, "wchar"))
        record.native_id = (cell(row, "native_id") or "").strip()
        record.queue = (cell(row, "queue") or "").strip()
        raw_action = cell(row, "error_action")
        record.error_action = "" if _is_null(raw_action) else raw_action.strip()

        if status == STATUS_COMPLETED:
            result.completed.append(record)
        else:
            result.failures.append(record)

    return result


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

class ProcessStats:
    """Everything observed for one process, pooled across every run."""

    def __init__(self, name: str, tier: str, module_process: Optional[str]) -> None:
        self.name = name
        self.tier = tier
        self.module_process = module_process
        self.tasks: List[TaskRecord] = []
        self.failures: List[TaskRecord] = []
        self.fits: Dict[str, Optional[FitResult]] = {}
        self.rec_cpus: Optional[int] = None
        self.rec_memory: Optional[float] = None
        self.rec_time: Optional[float] = None
        self.basis = "none"
        self.forced_by_failure = False
        self.clamped: List[str] = []

    # ── Observations ────────────────────────────────────────────────────────

    @property
    def n_tasks(self) -> int:
        return len(self.tasks)

    @property
    def run_keys(self) -> List[str]:
        return sorted({task.run_key for task in self.tasks})

    @property
    def n_runs(self) -> int:
        return len(self.run_keys)

    @property
    def n_retries(self) -> int:
        """Tasks that only succeeded on a second or later attempt."""
        return sum(1 for task in self.tasks if task.attempt > 1)

    @property
    def first_attempt_tasks(self) -> List[TaskRecord]:
        """Tasks whose request was not inflated by ``task.attempt``.

        Efficiency ratios use these only: a retry's request is a multiple of the
        declared value, so its ratio measures the retry ladder rather than the
        tool's appetite.
        """
        return [task for task in self.tasks if task.attempt <= 1]

    @property
    def max_peak_rss(self) -> Optional[float]:
        values = [task.peak_rss for task in self.tasks if task.peak_rss is not None]
        return max(values) if values else None

    @property
    def min_peak_rss(self) -> Optional[float]:
        """Smallest peak observed, which is the floor a dynamic request may fall to."""
        values = [task.peak_rss for task in self.tasks if task.peak_rss is not None]
        return min(values) if values else None

    @property
    def p95_realtime(self) -> Optional[float]:
        return percentile([task.realtime for task in self.tasks
                           if task.realtime is not None], 95)

    @property
    def max_realtime(self) -> Optional[float]:
        values = [task.realtime for task in self.tasks if task.realtime is not None]
        return max(values) if values else None

    @property
    def p90_cpu_cores(self) -> Optional[float]:
        return percentile([task.cpu_cores_used for task in self.tasks
                           if task.cpu_cores_used is not None], 90)

    @property
    def current_cpus(self) -> Optional[float]:
        values = [task.req_cpus for task in self.first_attempt_tasks
                  if task.req_cpus is not None]
        return median(values) if values else None

    @property
    def current_memory(self) -> Optional[float]:
        values = [task.req_memory for task in self.first_attempt_tasks
                  if task.req_memory is not None]
        return median(values) if values else None

    @property
    def current_time(self) -> Optional[float]:
        values = [task.req_time for task in self.first_attempt_tasks
                  if task.req_time is not None]
        return median(values) if values else None

    @property
    def memory_efficiency(self) -> Optional[float]:
        if self.max_peak_rss is None or not self.current_memory:
            return None
        return self.max_peak_rss / self.current_memory

    @property
    def cpu_efficiency(self) -> Optional[float]:
        values = [task.cpu_efficiency() for task in self.first_attempt_tasks]
        values = [value for value in values if value is not None]
        return sum(values) / len(values) if values else None

    @property
    def time_efficiency(self) -> Optional[float]:
        if self.p95_realtime is None or not self.current_time:
            return None
        return self.p95_realtime / self.current_time

    @property
    def cpu_hours(self) -> float:
        """Total allocated core-hours, used to rank what is worth tuning."""
        total = 0.0
        for task in self.tasks:
            if task.realtime is not None and task.req_cpus:
                total += task.realtime * task.req_cpus / 3600.0
        return total

    def size_points(self, metric: str) -> List[Tuple[float, float, str]]:
        """``(rchar, metric, run_key)`` triples for the scaling fit and scatters."""
        points: List[Tuple[float, float, str]] = []
        for task in self.tasks:
            if task.rchar is None or task.rchar <= 0:
                continue
            value = {
                "memory": task.peak_rss,
                "time": task.realtime,
                "cpu": task.cpu_cores_used,
            }.get(metric)
            if value is not None and value > 0:
                points.append((task.rchar, value, task.run_key))
        return points


def backfill_requests(tasks: Sequence[TaskRecord],
                      tiers: Dict[str, TierSpec]) -> List[str]:
    """Fill missing requested resources from the tier each task's label declares.

    Nextflow's *default* ``trace.fields`` records only what a task used -- no
    ``cpus``, ``memory`` or ``time``. Any trace written before this pipeline set an
    explicit field list, or with a bare ``-with-trace``, therefore has no
    denominator, and every efficiency figure would be blank.

    The declared value from ``conf/base.config`` is the best stand-in available, and
    for most processes the tier value *is* what was requested on the first attempt.
    It is an approximation in three ways, so it is recorded as inferred and reported
    as such rather than passed off as measured: it cannot see ``task.attempt``
    escalation, nor clamping by the site profile, and the handful of modules that
    carry an explicit ``memory { BcaResources.scaledMemory(...) }`` directive
    (``STARSOLO_ALIGN``, ``STARSOLO_INDEX``, ``MTX_TO_H5AD``, ``MTX_TO_10X``,
    ``SATURATION_TABLE``, ``SUBSET_VELOCYTO_MATRICES``, ``VELOCITY_H5AD``) override
    their label's memory entirely, so their backfilled figure is the tier's flat
    value rather than the size-scaled one they actually ran with.

    Parameters
    ----------
    tasks:
        Task records to fill in place.
    tiers:
        Declared values, from :func:`parse_base_config`.

    Returns
    -------
    Names of the fields that were backfilled for at least one task.
    """
    filled: set = set()
    for task in tasks:
        spec = tiers.get(task.tier or DEFAULT_TIER) or tiers.get(DEFAULT_TIER)
        if spec is None:
            continue
        if task.req_cpus is None and spec.cpus:
            task.req_cpus = spec.cpus
            filled.add("cpus")
        if task.req_memory is None and spec.memory:
            task.req_memory = spec.memory
            filled.add("memory")
        if task.req_time is None and spec.time:
            task.req_time = spec.time
            filled.add("time")
    return sorted(filled)


def aggregate(all_tasks: Sequence[TaskRecord],
              all_failures: Sequence[TaskRecord]) -> Dict[str, ProcessStats]:
    """Group task records by process name."""
    stats: Dict[str, ProcessStats] = {}
    for task in all_tasks:
        key = simple_process_name(task.process, task.name)
        if key not in stats:
            stats[key] = ProcessStats(key, task.tier or DEFAULT_TIER, task.module_process)
        stats[key].tasks.append(task)
    for task in all_failures:
        key = simple_process_name(task.process, task.name)
        if key not in stats:
            stats[key] = ProcessStats(key, task.tier or DEFAULT_TIER, task.module_process)
        stats[key].failures.append(task)
    return stats


# ---------------------------------------------------------------------------
# Recommendation
# ---------------------------------------------------------------------------

def _round_up(value: float, step: float) -> float:
    return math.ceil(value / step) * step


def _round_memory(value: float) -> float:
    """Quantise a memory recommendation, finely while it is small."""
    step = (MEMORY_FINE_STEP_BYTES if value < MEMORY_FINE_BELOW_BYTES
            else MEMORY_STEP_BYTES)
    return _round_up(value, step)


def _round_time(value: float) -> float:
    """Quantise a walltime recommendation, coarsely once it is long."""
    step = (TIME_COARSE_STEP_SECONDS if value > TIME_COARSE_ABOVE_SECONDS
            else TIME_STEP_SECONDS)
    return _round_up(value, step)


def recommend(stats: ProcessStats, args: argparse.Namespace) -> None:
    """Fill in ``rec_cpus``/``rec_memory``/``rec_time`` on *stats*.

    The three resources get deliberately different policies, because their
    failure modes differ:

    * **memory** -- a hard cliff (OOM kills the task), and per-process sample
      counts are small, so the maximum is used rather than a percentile.
    * **time** -- also a cliff, but over-requesting is nearly free, so p95 with a
      generous factor avoids one I/O-stalled task setting the tier.
    * **cpus** -- over-requesting only wastes allocation, never fails, so this is
      the one resource recommended downward, and it is never raised above what is
      already requested.

    When a trusted scaling fit exists and ``--target-size`` was given, the
    prediction at that size replaces the flat statistic.
    """
    if stats.n_tasks < args.min_tasks:
        stats.basis = f"insufficient data (n={stats.n_tasks})"
        stats.rec_cpus = int(stats.current_cpus) if stats.current_cpus else None
        stats.rec_memory = stats.current_memory
        stats.rec_time = stats.current_time
        return

    for metric in ("memory", "time", "cpu"):
        stats.fits[metric] = fit_loglog(stats.size_points(metric))

    memory_fit = stats.fits.get("memory")
    time_fit = stats.fits.get("time")
    use_fit = bool(args.target_size) and memory_fit is not None and memory_fit.trusted
    target_bytes = (args.target_size * 1024 ** 3) if args.target_size else None

    # ── Memory ──────────────────────────────────────────────────────────────
    observed_memory = stats.max_peak_rss
    if use_fit and target_bytes:
        predicted = memory_fit.predict(target_bytes)
        observed_memory = max(observed_memory or 0.0, predicted)
        stats.basis = (f"fit@{args.target_size:g}GB b={memory_fit.slope:.2f} "
                       f"r2={memory_fit.r2:.2f}")
    elif memory_fit is not None and memory_fit.trusted:
        # The fit is sound but nothing asked for an extrapolation, so the number
        # is still the observed maximum. Say so, otherwise the exponent column
        # looks like it fed a recommendation it did not.
        stats.basis = f"flat (b={memory_fit.slope:.2f} usable)"
    else:
        stats.basis = f"flat (n={stats.n_tasks}, runs={stats.n_runs})"

    if observed_memory:
        recommended = max(observed_memory * args.safety_memory, MIN_MEMORY_BYTES)
        stats.rec_memory = _round_memory(recommended)
    else:
        stats.rec_memory = stats.current_memory

    # ── Time ────────────────────────────────────────────────────────────────
    observed_time = stats.p95_realtime
    if args.target_size and time_fit is not None and time_fit.trusted and target_bytes:
        observed_time = max(observed_time or 0.0, time_fit.predict(target_bytes))
    if observed_time:
        recommended = max(observed_time * args.safety_time, MIN_TIME_SECONDS)
        stats.rec_time = _round_time(recommended)
    else:
        stats.rec_time = stats.current_time

    # ── CPUs ────────────────────────────────────────────────────────────────
    cores = stats.p90_cpu_cores
    if cores is not None:
        proposed = max(1, math.ceil(cores))
        if stats.current_cpus:
            proposed = min(proposed, int(stats.current_cpus))
        stats.rec_cpus = proposed
    elif stats.current_cpus:
        stats.rec_cpus = int(stats.current_cpus)

    # ── Observed kills override every statistic ─────────────────────────────
    # A task the scheduler killed proves its request was too small, whatever the
    # surviving tasks suggest, so the recommendation must clear it.
    #
    # How far to clear it depends on what else is known. When a later attempt
    # succeeded, its peak_rss is a real measurement of what the task needs and the
    # percentile logic above has already used it -- the kill then only supplies a
    # floor. With no successful observation there is no measurement at all, and
    # doubling the request that died is the only defensible guess.
    killed = [task for task in stats.failures if task.killed]
    if killed:
        stats.forced_by_failure = True
        have_measurement = stats.max_peak_rss is not None
        factor = 1.1 if have_measurement else 2.0

        memory_requests = [task.req_memory for task in killed if task.req_memory]
        if memory_requests:
            stats.rec_memory = _round_memory(
                max(stats.rec_memory or 0.0, max(memory_requests) * factor))

        time_factor = 1.1 if stats.max_realtime is not None else 2.0
        time_requests = [task.req_time for task in killed if task.req_time]
        if time_requests:
            stats.rec_time = _round_time(
                max(stats.rec_time or 0.0, max(time_requests) * time_factor))

    # ── Never recommend below what was already observed ─────────────────────
    if stats.max_peak_rss and stats.rec_memory:
        stats.rec_memory = max(stats.rec_memory,
                               _round_memory(stats.max_peak_rss * 1.1))
    if stats.max_realtime and stats.rec_time:
        stats.rec_time = max(stats.rec_time, _round_time(stats.max_realtime * 1.1))

    # ── Site limits ─────────────────────────────────────────────────────────
    if args.max_memory and stats.rec_memory and stats.rec_memory > args.max_memory * 1024 ** 3:
        stats.rec_memory = args.max_memory * 1024 ** 3
        stats.clamped.append("memory")
    if args.max_time and stats.rec_time and stats.rec_time > args.max_time * 3600:
        stats.rec_time = args.max_time * 3600
        stats.clamped.append("time")
    if args.max_cpus and stats.rec_cpus and stats.rec_cpus > args.max_cpus:
        stats.rec_cpus = args.max_cpus
        stats.clamped.append("cpus")


# ---------------------------------------------------------------------------
# Tool-internal memory caps
# ---------------------------------------------------------------------------

class CapFinding:
    """A tool-internal memory ceiling that blocks or endangers a recommendation."""

    def __init__(self, process: str, param: str, note: str) -> None:
        self.process = process
        self.param = param
        self.note = note
        self.pinned_values: Dict[str, Optional[float]] = {}   # run key -> bytes
        self.severity = "info"
        self.message = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "process": self.process,
            "param": self.param,
            "note": self.note,
            "pinned_values": self.pinned_values,
            "severity": self.severity,
            "message": self.message,
        }


def check_tool_caps(stats_by_name: Dict[str, ProcessStats],
                    runs: Sequence[RunContext]) -> List[CapFinding]:
    """Flag processes whose recommendation a tool-internal cap would neutralise.

    A ``memory`` directive only sets what the scheduler reserves. STAR additionally
    takes its own ceilings on the command line, and those do not move when the
    Nextflow directive does. Three situations matter, in increasing severity:

    * the cap is left on auto -- fine, it now tracks ``task.memory``;
    * the cap is pinned **below** the recommendation -- the extra memory would be
      reserved and then refused by the tool, so the recommendation cannot take
      effect until the pin is raised or removed;
    * the cap is pinned **above** the recommendation -- the tool will try to use
      more than the job reserves and the cgroup will kill it.

    Parameters
    ----------
    stats_by_name:
        Per-process statistics, already carrying recommendations.
    runs:
        Runs whose ``run_config_*.txt`` supplies the recorded parameter values.

    Returns
    -------
    One :class:`CapFinding` per affected process that was actually observed.
    """
    findings: List[CapFinding] = []

    for process, (param, note) in TOOL_MEMORY_CAPS.items():
        stats = stats_by_name.get(process)
        if stats is None or not stats.n_tasks:
            continue

        finding = CapFinding(process, param, note)
        pinned: List[float] = []
        for run in runs:
            raw = run.params.get(param)
            if raw is None:
                finding.pinned_values[run.key] = None
                continue
            token = str(raw).strip().strip('"').strip("'").lower()
            if token in CAP_AUTO_TOKENS:
                finding.pinned_values[run.key] = None
                continue
            value = parse_bytes(token)
            finding.pinned_values[run.key] = value
            if value:
                pinned.append(value)

        recommended = stats.rec_memory
        if not pinned:
            finding.severity = "ok"
            finding.message = (
                f"{param} was left on auto in every run, so it is derived from "
                f"task.memory in the module and follows any change to the "
                f"{stats.tier} tier. Nothing to do.")
            findings.append(finding)
            continue

        lowest, highest = min(pinned), max(pinned)
        if recommended and highest > recommended:
            finding.severity = "danger"
            finding.message = (
                f"{param} is pinned at {fmt_bytes(highest)}, above the recommended "
                f"request of {fmt_bytes(recommended)}. The tool would try to use more "
                f"than the job reserves and the cgroup would kill it. Either raise the "
                f"request to at least {fmt_bytes(highest)} or unset {param} so it is "
                f"derived from task.memory.")
        elif recommended and lowest < recommended * 0.75:
            finding.severity = "warn"
            finding.message = (
                f"{param} is pinned at {fmt_bytes(lowest)}, well below the recommended "
                f"request of {fmt_bytes(recommended)}. Raising the request alone will "
                f"not help -- the tool will still refuse to use more than its own "
                f"ceiling. Unset {param} so it is derived from task.memory, or raise it "
                f"alongside the request.")
        else:
            finding.severity = "info"
            finding.message = (
                f"{param} is pinned at {fmt_bytes(lowest)}, roughly in line with the "
                f"recommended request. Remember it does not move when the tier does.")
        findings.append(finding)

    return findings


# ---------------------------------------------------------------------------
# Config emission
# ---------------------------------------------------------------------------

class TierPlan:
    """A retuned ``withLabel:`` block plus the members lifted out of it."""

    def __init__(self, label: str) -> None:
        self.label = label
        self.members: List[ProcessStats] = []
        self.outliers: List[ProcessStats] = []
        self.cpus: Optional[int] = None
        self.memory: Optional[float] = None
        self.time: Optional[float] = None
        self.changed = False
        self.n_declared = 0        # processes in modules/ carrying this label
        self.capped_by_coverage: List[str] = []


def build_tier_plans(stats_by_name: Dict[str, ProcessStats],
                     tiers: Dict[str, TierSpec],
                     module_labels: Dict[str, ProcessLabels],
                     args: argparse.Namespace) -> Dict[str, TierPlan]:
    """Split each tier into a shared recommendation plus per-process outliers.

    ``process_medium`` covers 19 modules here, so sizing the tier to its hungriest
    member would over-allocate the other 18 on every task. Any member needing more
    than :data:`OUTLIER_FACTOR` times the tier median is therefore emitted as its
    own ``withName:`` block, and the tier is sized from what remains.

    A tier is only ever *lowered* when every process declaring it was actually
    observed. A run exercises whichever branch its protocol selects, so a typical
    trace covers a handful of a tier's members; shrinking the shared value on that
    evidence would silently under-provision every member that did not run, and the
    first sign of it would be an OOM in an unrelated pipeline branch. Raising a
    tier is always safe, so that is allowed on partial coverage.
    """
    plans: Dict[str, TierPlan] = {}
    declared_counts: Dict[str, int] = {}
    for info in module_labels.values():
        declared_counts[info.tier] = declared_counts.get(info.tier, 0) + 1

    for label in list(RESOURCE_LABELS) + [DEFAULT_TIER]:
        plan = TierPlan(label)
        plan.n_declared = declared_counts.get(label, 0)
        eligible = [stats for stats in stats_by_name.values()
                    if stats.tier == label and stats.n_tasks >= args.min_tasks]
        if not eligible:
            plans[label] = plan
            continue

        # The unlabelled default block is shared with anything that never declared
        # a label, including future modules, so it is never widened as a tier --
        # its members are always emitted individually.
        if label == DEFAULT_TIER:
            plan.outliers = eligible
            plans[label] = plan
            continue

        memory_values = [stats.rec_memory for stats in eligible if stats.rec_memory]
        time_values = [stats.rec_time for stats in eligible if stats.rec_time]
        median_memory = median(memory_values) if memory_values else None
        median_time = median(time_values) if time_values else None

        outliers: List[ProcessStats] = []
        for stats in eligible:
            # A killed task is a per-process fact, not a statistical one, so it
            # lifts the process out regardless of how many members the tier has:
            # its raised request must not be imposed on every sibling.
            if stats.forced_by_failure:
                outliers.append(stats)
                continue
            if len(eligible) < MIN_MEMBERS_FOR_OUTLIER:
                continue
            over_memory = (median_memory and stats.rec_memory
                           and stats.rec_memory > median_memory * OUTLIER_FACTOR)
            over_time = (median_time and stats.rec_time
                         and stats.rec_time > median_time * OUTLIER_FACTOR)
            if over_memory or over_time:
                outliers.append(stats)

        members = [stats for stats in eligible if stats not in outliers]
        plan.members = members
        plan.outliers = outliers

        pool = members or eligible
        cpu_values = [stats.rec_cpus for stats in pool if stats.rec_cpus]
        memory_pool = [stats.rec_memory for stats in pool if stats.rec_memory]
        time_pool = [stats.rec_time for stats in pool if stats.rec_time]
        plan.cpus = max(cpu_values) if cpu_values else None
        plan.memory = max(memory_pool) if memory_pool else None
        plan.time = max(time_pool) if time_pool else None

        current = tiers.get(label)

        # Partial coverage: never shrink a shared value on the strength of the few
        # members this trace happened to exercise.
        #
        # Count in the same namespace on both sides. n_declared counts *module*
        # process names, so the observed set has to be folded back to those too:
        # an aliased inclusion appears in the trace under each alias
        # (DOUBLET_FILTER_RAW, DOUBLET_FILTER_CELL_CALLED) but is one declared
        # module (DOUBLET_FILTER). Counting the aliases would make unobserved
        # negative and silently switch this guard off for the whole tier.
        observed_names = {stats.module_process or stats.name for stats in eligible}
        unobserved = plan.n_declared - len(observed_names)
        if current and unobserved > 0:
            if plan.cpus and current.cpus and plan.cpus < current.cpus:
                plan.cpus = int(current.cpus)
                plan.capped_by_coverage.append("cpus")
            if plan.memory and current.memory and plan.memory < current.memory:
                plan.memory = current.memory
                plan.capped_by_coverage.append("memory")
            if plan.time and current.time and plan.time < current.time:
                plan.time = current.time
                plan.capped_by_coverage.append("time")
        if current:
            def near(new: Optional[float], old: Optional[float]) -> bool:
                if new is None or not old:
                    return True
                return abs(new - old) / old <= NO_CHANGE_TOLERANCE

            plan.changed = not (near(plan.cpus, current.cpus)
                                and near(plan.memory, current.memory)
                                and near(plan.time, current.time))
        else:
            plan.changed = True

        plans[label] = plan

    return plans


def render_config(plans: Dict[str, TierPlan],
                  tiers: Dict[str, TierSpec],
                  runs: Sequence[RunContext],
                  n_tasks: int,
                  cap_findings: Sequence[CapFinding],
                  args: argparse.Namespace) -> str:
    """Render the tuned ``withLabel``/``withName`` blocks as a Nextflow config.

    Every block carries its provenance -- task and run counts, the observed
    maximum, the scaling exponent and whether the basis was a fit or a flat
    statistic -- because a resource number a reader cannot audit will not be
    trusted, and rightly so.
    """
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines: List[str] = [
        "/*",
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
        "    Tuned resource requests -- GENERATED FILE, DO NOT EDIT BY HAND",
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
        f"    Generated by bin/resource_efficiency.py on {stamp}",
        f"    Derived from {len(runs)} run(s), {n_tasks} completed task(s):",
    ]
    for run in runs:
        lines.append(f"      {run.trace_path}")
    lines += [
        "",
        "    Overrides conf/base.config when passed after it. submit_nextflow.sh adds",
        "    it automatically when present; otherwise:",
        "",
        "        nextflow run ... -c conf/resources_tuned.config",
        "",
        "    Regenerate after a few more runs; delete the file to revert to base.config.",
    ]

    actionable = [finding for finding in cap_findings if finding.severity in ("warn", "danger")]
    if actionable:
        lines += [
            "",
            "    !! These requests are NOT sufficient on their own. The tools below take",
            "    !! their own memory ceilings, which do not move when this file changes:",
        ]
        for finding in actionable:
            lines.append(f"    !!   {finding.process}: {finding.param} -- {finding.message}")

    lines += [
        "----------------------------------------------------------------------------------------",
        "*/",
        "",
        "process {",
    ]

    def emit_directives(indent: str, cpus: Optional[int], memory: Optional[float],
                        time_value: Optional[float]) -> List[str]:
        # One padding width for all three, so the block lines up the way
        # conf/base.config does and a diff between them stays readable.
        width = 7
        out: List[str] = []
        if cpus:
            out.append(f"{indent}cpus   = {{ {str(cpus):<{width}} * task.attempt }}")
        if memory:
            out.append(f"{indent}memory = {{ {groovy_memory(memory):<{width}} * task.attempt }}")
        if time_value:
            out.append(f"{indent}time   = {{ {groovy_time(time_value):<{width}} * task.attempt }}")
        return out

    emitted = False
    for label in RESOURCE_LABELS:
        plan = plans.get(label)
        if not plan or not plan.members:
            continue
        # A tier the coverage guard pinned all the way back to its conf/base.config
        # value comes out 'unchanged', which used to drop the whole block -- hiding
        # the guard precisely when it did the most work, and leaving the reader of
        # this config no way to tell a tier was deliberately held from one that was
        # never looked at. Such a tier is reported below as comments, with no
        # withLabel block, since there is nothing to override.
        if not plan.changed and not plan.capped_by_coverage:
            continue
        emitted = emitted or plan.changed
        current = tiers.get(label)
        n_member_tasks = sum(stats.n_tasks for stats in plan.members)
        n_member_runs = len({key for stats in plan.members for key in stats.run_keys})
        lines.append("")
        lines.append(f"    // {label}: {len(plan.members)} process(es), "
                     f"{n_member_tasks} task(s), {n_member_runs} run(s).")
        if current:
            lines.append(f"    // was: cpus={current.cpus and int(current.cpus)} "
                         f"memory={fmt_bytes(current.memory)} time={fmt_duration(current.time)}")
        worst = max(plan.members, key=lambda s: s.max_peak_rss or 0.0)
        lines.append(f"    // driven by {worst.name}: max peak_rss {fmt_bytes(worst.max_peak_rss)}"
                     f" ({fmt_pct(worst.memory_efficiency)} of request); {worst.basis}")
        if plan.outliers:
            lines.append(f"    // lifted out into withName blocks below: "
                         f"{', '.join(stats.name for stats in plan.outliers)}")
        # Same namespace fold as in plan_tiers(), so the reported count matches the
        # one the coverage guard actually acted on.
        unobserved = plan.n_declared - len(
            {stats.module_process or stats.name for stats in plan.members}
            | {stats.module_process or stats.name for stats in plan.outliers})
        if unobserved > 0:
            lines.append(f"    // {unobserved} of {plan.n_declared} process(es) with this label "
                         f"were not exercised by these runs.")
            if plan.capped_by_coverage:
                lines.append(f"    // {'/'.join(plan.capped_by_coverage)} therefore held at the "
                             f"conf/base.config value rather than lowered, so the")
                lines.append("    // processes that never ran are not silently under-provisioned.")
        if plan.changed:
            lines.append(f"    withLabel:{label} {{")
            lines += emit_directives("        ", plan.cpus, plan.memory, plan.time)
            lines.append("    }")
        else:
            lines.append("    // Nothing to override: every dimension ended up back at the "
                         "conf/base.config value.")

    for label in list(RESOURCE_LABELS) + [DEFAULT_TIER]:
        plan = plans.get(label)
        if not plan:
            continue
        for stats in sorted(plan.outliers, key=lambda s: -(s.max_peak_rss or 0.0)):
            emitted = True
            origin = ("unlabelled, inherits the shared process{} default"
                      if label == DEFAULT_TIER else f"outlier within {label}")
            lines.append("")
            lines.append(f"    // {stats.name}: {origin}. "
                         f"{stats.n_tasks} task(s), {stats.n_runs} run(s).")
            lines.append(f"    // max peak_rss {fmt_bytes(stats.max_peak_rss)}"
                         f" of {fmt_bytes(stats.current_memory)} requested"
                         f" ({fmt_pct(stats.memory_efficiency)}); {stats.basis}")
            if stats.n_retries:
                lines.append(f"    // {stats.n_retries} task(s) needed a retry.")
            if stats.forced_by_failure:
                lines.append("    // raised: the scheduler killed at least one attempt.")
            if stats.clamped:
                lines.append(f"    // clamped to --max-{'/'.join(stats.clamped)}.")
            lines.append(f"    withName:'.*:{stats.name}' {{")
            lines += emit_directives("        ", stats.rec_cpus, stats.rec_memory, stats.rec_time)
            lines.append("    }")

    if not emitted:
        lines.append("")
        lines.append("    // Nothing to change: every tier is already within "
                     f"{int(NO_CHANGE_TOLERANCE * 100)}% of its observed requirement,")
        lines.append("    // or too few tasks were observed to justify a change.")

    lines.append("}")
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Dynamic (input-size dependent) directives
# ---------------------------------------------------------------------------

def render_dynamic(ordered: Sequence[ProcessStats],
                   tiers: Dict[str, TierSpec],
                   args: argparse.Namespace) -> str:
    """Render a ``params.dynamic_memory`` block sizing processes from their input.

    A flat tier has to cover the largest dataset it will ever see, so on everything
    smaller the difference is reserved and left idle. Where usage provably tracks
    input volume, ``lib/BcaResources.groovy`` scales the request with it instead.

    The fit gives ``memory = 10**intercept * bytes**exponent``. That is re-expressed
    here as an anchor -- "at ``ref_gb`` of input the process needed ``mem_gb``" --
    because an anchor is something a reader can sanity-check against their own data,
    where a raw coefficient is not. ``ref_gb`` is the geometric mean of the observed
    input sizes, which is the centre of the fitted range and so the point the fit is
    most confident about.

    Only a process whose fit passed every quality gate gets an entry; the rest keep
    a flat request, which is what the evidence supports.
    """
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines: List[str] = [
        "/*",
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
        "    Input-size dependent memory -- GENERATED FILE, DO NOT EDIT BY HAND",
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
        f"    Generated by bin/resource_efficiency.py on {stamp}",
        "",
        "    Replaces the `dynamic_memory` map in nextflow.config. A process only",
        "    honours its entry if its module calls BcaResources.scaledMemory(); the",
        "    entries below for processes that do not are inert until one does.",
        "",
        "    Reading: at `ref_gb` of input the process needed `mem_gb`, and memory",
        "    grows as size**`exponent` from there, clamped to [floor_gb, cap_gb].",
        "",
        "    Apply with:  nextflow run ... -c <this file>",
        "----------------------------------------------------------------------------------------",
        "*/",
        "",
        "params {",
        "    dynamic_memory = [",
    ]

    emitted = 0
    for stats in ordered:
        fit = stats.fits.get("memory")
        if fit is None or not fit.trusted or not stats.rec_memory:
            continue

        sizes = [task.rchar for task in stats.tasks if task.rchar and task.rchar > 0]
        if not sizes:
            continue
        emitted += 1

        # Anchor at the geometric mean of the observed sizes: on log axes that is
        # the centre of the fitted range, so it is where the fit is best supported.
        log_mean = sum(math.log10(size) for size in sizes) / len(sizes)
        ref_bytes = 10 ** log_mean
        anchor_bytes = (10 ** (fit.intercept + fit.slope * log_mean)) * args.safety_memory

        ref_gb = max(0.001, ref_bytes / 1024 ** 3)
        mem_gb = max(1, int(math.ceil(anchor_bytes / 1024 ** 3)))

        # The floor is the *smallest* peak ever observed, not a fraction of the
        # largest: a floor near the maximum would hand every small input the flat
        # request back and undo the point of sizing dynamically.
        floor_gb = max(1, int(math.ceil(
            (stats.min_peak_rss or 0.0) * args.safety_memory / 1024 ** 3)))

        # The cap has to clear both the largest prediction and the declared tier --
        # a process already close to exhausting its tier needs room above it, or the
        # cap simply reinstates the failure the scaling was meant to avoid.
        predicted_hi = (10 ** (fit.intercept + fit.slope * math.log10(max(sizes)))) \
            * args.safety_memory
        tier = tiers.get(stats.tier)
        tier_gb = (tier.memory / 1024 ** 3) if tier and tier.memory else 0
        cap_gb = max(1, int(math.ceil(max(predicted_hi * 2 / 1024 ** 3, tier_gb))))

        note = ""
        if fit.slope < 0.15:
            note = "  // exponent near zero: barely scales, a flat request is nearly as good"

        lines += [
            f"        // {stats.name}: {stats.n_tasks} task(s), {stats.n_runs} run(s); "
            f"exponent {fit.slope:.2f}, r2 {fit.r2:.2f}.",
            f"        //   observed input {fmt_bytes(min(sizes))} to {fmt_bytes(max(sizes))}; "
            f"peak memory {fmt_bytes(stats.max_peak_rss)}; flat would be "
            f"{fmt_bytes(stats.rec_memory)}.{note}",
            f"        {stats.name}: [ exponent: {fit.slope:.2f}, ref_gb: {ref_gb:.3g}, "
            f"mem_gb: {mem_gb}, floor_gb: {floor_gb}, cap_gb: {cap_gb} ],",
        ]

    if not emitted:
        lines.append("        // Nothing qualified: no process had a scaling fit that passed")
        lines.append("        // every quality gate, so sizing from input would be guesswork.")

    lines += ["    ]", "}", ""]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Text report
# ---------------------------------------------------------------------------

def render_table(ordered: Sequence[ProcessStats]) -> str:
    """Render the per-process efficiency and recommendation table."""
    header = (f"{'process':<32} {'label':<20} {'runs':>4} {'task':>5} {'retr':>4} "
              f"{'memory used/req':>22} {'eff':>5} {'cpu':>5} {'time used/req':>20} "
              f"{'b':>6} {'-> recommendation':<34} basis")
    lines = [header, "-" * len(header)]

    for stats in ordered:
        memory_pair = f"{fmt_bytes(stats.max_peak_rss)}/{fmt_bytes(stats.current_memory)}"
        time_pair = f"{fmt_duration(stats.p95_realtime)}/{fmt_duration(stats.current_time)}"
        fit = stats.fits.get("memory")
        exponent = f"{fit.slope:.2f}" if fit and fit.trusted else "-"
        recommendation = " ".join(filter(None, [
            f"{stats.rec_cpus}cpu" if stats.rec_cpus else "",
            fmt_bytes(stats.rec_memory) if stats.rec_memory else "",
            fmt_duration(stats.rec_time) if stats.rec_time else "",
        ]))
        flags = ""
        if stats.forced_by_failure:
            flags += " [killed]"
        if stats.clamped:
            flags += " [clamped]"
        lines.append(
            f"{stats.name[:32]:<32} {stats.tier[:20]:<20} {stats.n_runs:>4} "
            f"{stats.n_tasks:>5} {stats.n_retries:>4} {memory_pair:>22} "
            f"{fmt_pct(stats.memory_efficiency):>5} {fmt_pct(stats.cpu_efficiency):>5} "
            f"{time_pair:>20} {exponent:>6} {recommendation:<34} {stats.basis}{flags}")

    return "\n".join(lines)


def write_tsv(path: str, ordered: Sequence[ProcessStats]) -> None:
    """Write the machine-readable per-process summary."""
    columns = [
        "process", "label", "n_runs", "n_tasks", "n_retries", "n_failures",
        "cpu_hours", "req_cpus", "req_memory_bytes", "req_time_seconds",
        "max_peak_rss_bytes", "p95_realtime_seconds", "p90_cpu_cores",
        "memory_efficiency", "cpu_efficiency", "time_efficiency",
        "memory_exponent", "memory_fit_intercept", "memory_fit_r2", "memory_fit_trusted",
        "rec_cpus", "rec_memory_bytes", "rec_time_seconds", "basis",
        "forced_by_failure", "clamped",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        for stats in ordered:
            fit = stats.fits.get("memory")
            writer.writerow([
                stats.name, stats.tier, stats.n_runs, stats.n_tasks,
                stats.n_retries, len(stats.failures), f"{stats.cpu_hours:.2f}",
                stats.current_cpus, stats.current_memory, stats.current_time,
                stats.max_peak_rss, stats.p95_realtime, stats.p90_cpu_cores,
                stats.memory_efficiency, stats.cpu_efficiency, stats.time_efficiency,
                fit.slope if fit else "", fit.intercept if fit else "",
                fit.r2 if fit else "", fit.trusted if fit else "",
                stats.rec_cpus, stats.rec_memory, stats.rec_time, stats.basis,
                stats.forced_by_failure, ",".join(stats.clamped),
            ])


def build_json(ordered: Sequence[ProcessStats],
               runs: Sequence[RunContext],
               unmapped: Sequence[str],
               parse_results: Sequence[TraceParseResult],
               cap_findings: Sequence[CapFinding]) -> Dict[str, Any]:
    """Assemble the machine-readable dump used by the test check."""
    return {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "runs": [
            {"key": run.key, "trace": run.trace_path, "results_dir": run.results_dir,
             "n_samples": len(run.samples),
             "protocol": run.params.get("protocol"),
             "mapping_software": run.params.get("mapping_software")}
            for run in runs
        ],
        "n_completed_tasks": sum(len(result.completed) for result in parse_results),
        "n_failed_tasks": sum(len(result.failures) for result in parse_results),
        "n_cached_rows": sum(result.n_cached for result in parse_results),
        "raw_mode": any(result.raw_mode for result in parse_results),
        "missing_fields": sorted({field for result in parse_results
                                  for field in result.missing_fields}),
        "unmapped_processes": list(unmapped),
        "tool_memory_caps": [finding.to_dict() for finding in cap_findings],
        "processes": [
            {
                "process": stats.name,
                "label": stats.tier,
                "module_process": stats.module_process,
                "n_runs": stats.n_runs,
                "n_tasks": stats.n_tasks,
                "n_retries": stats.n_retries,
                "n_failures": len(stats.failures),
                "cpu_hours": stats.cpu_hours,
                "current_cpus": stats.current_cpus,
                "current_memory": stats.current_memory,
                "current_time": stats.current_time,
                "max_peak_rss": stats.max_peak_rss,
                "p95_realtime": stats.p95_realtime,
                "p90_cpu_cores": stats.p90_cpu_cores,
                "memory_efficiency": stats.memory_efficiency,
                "cpu_efficiency": stats.cpu_efficiency,
                "time_efficiency": stats.time_efficiency,
                "fits": {name: (fit.to_dict() if fit else None)
                         for name, fit in stats.fits.items()},
                "recommended_cpus": stats.rec_cpus,
                "recommended_memory": stats.rec_memory,
                "recommended_time": stats.rec_time,
                "basis": stats.basis,
                "forced_by_failure": stats.forced_by_failure,
                "clamped": stats.clamped,
            }
            for stats in ordered
        ],
        "failures": [
            {"process": simple_process_name(task.process, task.name),
             "run": task.run_key, "tag": task.tag, "exit": task.exit_code,
             "attempt": task.attempt, "req_memory": task.req_memory,
             "req_time": task.req_time, "killed": task.killed,
             "native_id": task.native_id}
            for stats in ordered for task in stats.failures
        ],
    }


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def encode_image(image_path: Optional[str]) -> Optional[str]:
    """Return a PNG data-URI for HTML embedding, or ``None`` if unreadable."""
    if not image_path or not os.path.exists(image_path):
        return None
    try:
        with open(image_path, "rb") as handle:
            encoded = base64.b64encode(handle.read()).decode("utf-8")
        return f"data:image/png;base64,{encoded}"
    except OSError as exc:
        sys.stderr.write(f"Warning: could not encode image {image_path}: {exc}\n")
        return None


def make_plots(ordered: Sequence[ProcessStats],
               runs: Sequence[RunContext],
               figures_dir: str,
               args: argparse.Namespace) -> Dict[str, str]:
    """Draw the overview charts and the per-process scatter panels.

    Returns
    -------
    ``{figure_key: png_path}``. Empty when matplotlib is unavailable, which is
    reported but never fatal -- the table and config emission do not need it.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")  # headless: this normally runs on a login node
        import matplotlib.pyplot as plt
    except ImportError:
        sys.stderr.write("Warning: matplotlib unavailable, skipping plots. "
                         "Pass --no-plots to silence this.\n")
        return {}

    os.makedirs(figures_dir, exist_ok=True)
    figures: Dict[str, str] = {}
    colours = plt.get_cmap("tab10")
    run_colour = {run.key: colours(i % 10) for i, run in enumerate(runs)}

    # ── Overview A: efficiency by process ───────────────────────────────────
    # Not capped by --max-panels: that limits the expensive per-process figures,
    # while this single chart is the one place every process stays visible.
    ranked = [stats for stats in ordered if stats.n_tasks]
    # Efficiency needs a requested value to divide by. When the trace recorded none
    # and no tier could be substituted, fall back to absolute usage: "STARSOLO_ALIGN
    # peaked at 96 GB" is still worth knowing, and an empty chart is not.
    have_efficiency = any(stats.memory_efficiency is not None
                          or stats.cpu_efficiency is not None for stats in ranked)
    labels = [f"{stats.name}{' *' if (stats.n_retries or stats.forced_by_failure) else ''}"
              for stats in ranked]

    if ranked and have_efficiency:
        height = max(3.0, 0.32 * len(ranked) + 2.0)
        fig, axis = plt.subplots(figsize=(11, height))
        positions = range(len(ranked))
        memory_eff = [(stats.memory_efficiency or 0) * 100 for stats in ranked]
        cpu_eff = [(stats.cpu_efficiency or 0) * 100 for stats in ranked]
        axis.barh([p + 0.2 for p in positions], memory_eff, height=0.4,
                  label="memory (max peak_rss / request)", color="#4C72B0")
        axis.barh([p - 0.2 for p in positions], cpu_eff, height=0.4,
                  label="cpu (mean cores used / request)", color="#DD8452")
        axis.set_yticks(list(positions))
        axis.set_yticklabels(labels, fontsize=8)
        axis.axvline(100, color="grey", linestyle="--", linewidth=1)
        axis.set_xlabel("% of requested resource actually used")
        axis.set_title("Resource efficiency by process  (* = had a retry or a kill)")
        axis.legend(loc="lower right", fontsize=8)
        axis.invert_yaxis()
        fig.tight_layout()
        path = os.path.join(figures_dir, "overview_efficiency.png")
        fig.savefig(path, dpi=110)
        plt.close(fig)
        figures["overview_efficiency"] = path

    elif ranked:
        # Absolute mode. Memory and cores need separate axes -- GB and core counts
        # share no scale, and putting them on one would make both unreadable.
        height = max(3.0, 0.32 * len(ranked) + 2.0)
        fig, axes = plt.subplots(1, 2, figsize=(13, height), sharey=True)
        positions = list(range(len(ranked)))
        axes[0].barh(positions, [(stats.max_peak_rss or 0) / 1024 ** 3 for stats in ranked],
                     height=0.6, color="#4C72B0")
        axes[0].set_xlabel("peak memory used (GB)")
        axes[0].set_yticks(positions)
        axes[0].set_yticklabels(labels, fontsize=8)
        axes[1].barh(positions, [stats.p90_cpu_cores or 0 for stats in ranked],
                     height=0.6, color="#DD8452")
        axes[1].set_xlabel("cores actually used (p90)")
        for axis in axes:
            axis.invert_yaxis()
            axis.grid(True, axis="x", alpha=0.15)
        fig.suptitle("Resource usage by process  (* = had a retry or a kill)\n"
                     "no requested values in this trace, so absolute usage is shown",
                     fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        path = os.path.join(figures_dir, "overview_efficiency.png")
        fig.savefig(path, dpi=110)
        plt.close(fig)
        figures["overview_efficiency"] = path

    # ── Overview B: label roll-up ───────────────────────────────────────────
    tier_members: Dict[str, List[ProcessStats]] = {}
    for stats in ordered:
        tier_members.setdefault(stats.tier, []).append(stats)
    if tier_members:
        fig, axis = plt.subplots(figsize=(11, max(3.0, 0.8 * len(tier_members) + 1.5)))
        tier_names = [label for label in list(RESOURCE_LABELS) + [DEFAULT_TIER]
                      if label in tier_members]
        drawn = 0
        for row, label in enumerate(tier_names):
            for stats in tier_members[label]:
                for task in stats.tasks:
                    # As with overview A, fall back to absolute peak_rss when there
                    # is no request to express it as a fraction of.
                    value = (task.memory_efficiency() * 100 if have_efficiency
                             and task.memory_efficiency() is not None
                             else (task.peak_rss / 1024 ** 3
                                   if task.peak_rss is not None else None))
                    if value is None:
                        continue
                    drawn += 1
                    axis.plot(value, row + (hash(stats.name) % 100 - 50) / 400.0,
                              "o", markersize=3, alpha=0.45,
                              color=run_colour.get(task.run_key, "#333333"))
        axis.set_yticks(range(len(tier_names)))
        axis.set_yticklabels(tier_names, fontsize=9)
        if have_efficiency:
            axis.axvline(100, color="crimson", linestyle="--", linewidth=1,
                         label="request (100%)")
            axis.set_xlabel("memory used as % of request, one point per task")
            axis.legend(loc="lower right", fontsize=8)
        else:
            axis.set_xlabel("peak memory used (GB), one point per task")
        axis.set_title("Per-tier spread: which processes drive each shared label")
        axis.invert_yaxis()
        fig.tight_layout()
        if drawn:
            path = os.path.join(figures_dir, "overview_tiers.png")
            fig.savefig(path, dpi=110)
            figures["overview_tiers"] = path
        else:
            # No memory was recorded at all; an axis with nothing on it tells the
            # reader less than its absence does.
            sys.stderr.write("Warning: no peak_rss values in the trace, "
                             "skipping the per-tier figure.\n")
        plt.close(fig)

    # ── Per-process scatter panels ──────────────────────────────────────────
    panels = [stats for stats in ordered if stats.n_tasks >= 2][:args.max_panels]
    for stats in panels:
        # Cores get a linear y axis: the values span a narrow band around the
        # request, and a log axis there turns "11 to 16 cores" into a misleading
        # near-flat line. Memory and runtime span decades and need the log axis.
        specs = [
            ("memory", stats.current_memory, stats.rec_memory,
             lambda v: v / 1024 ** 3, "peak RSS (GB)", "log"),
            ("time", stats.current_time, stats.rec_time,
             lambda v: v / 3600.0, "runtime (h)", "log"),
            ("cpu", stats.current_cpus, stats.rec_cpus,
             lambda v: v, "cores used (%cpu / 100)", "linear"),
        ]
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
        for axis, (metric, current, recommended, scale, ylabel, yscale) in zip(axes, specs):
            points = stats.size_points(metric)
            for run in runs:
                subset = [(x, y) for x, y, key in points if key == run.key]
                if not subset:
                    continue
                axis.scatter([x / 1024 ** 3 for x, _ in subset],
                             [scale(y) for _, y in subset],
                             s=26, alpha=0.75, label=run.key,
                             color=run_colour.get(run.key, "#333333"))

            fit = stats.fits.get(metric)
            if fit and fit.trusted and points:
                xs = sorted(x for x, _, _ in points)
                grid = [xs[0] * (xs[-1] / xs[0]) ** (i / 49.0) for i in range(50)] \
                    if xs[0] > 0 and xs[-1] > xs[0] else xs
                predicted = [10 ** (fit.intercept + fit.slope * math.log10(x)) for x in grid]
                upper = [value * 10 ** (2 * fit.resid_std) for value in predicted]
                lower = [value / 10 ** (2 * fit.resid_std) for value in predicted]
                gx = [x / 1024 ** 3 for x in grid]
                axis.plot(gx, [scale(v) for v in predicted], "-", color="#333333",
                          linewidth=1.2,
                          label=f"fit b={fit.slope:.2f} (r²={fit.r2:.2f})")
                axis.fill_between(gx, [scale(v) for v in lower], [scale(v) for v in upper],
                                  color="#333333", alpha=0.10)

            if current:
                axis.axhline(scale(current), color="grey", linestyle="--", linewidth=1.1,
                             label="current request")
            if recommended:
                axis.axhline(scale(recommended), color="#2CA02C", linestyle="--",
                             linewidth=1.1, label="recommended")

            # Killed tasks are plotted at the request they died under: the point
            # where the horizontal "current request" line was not enough.
            for task in stats.failures:
                if not task.killed or task.rchar is None:
                    continue
                value = {"memory": task.req_memory, "time": task.req_time,
                         "cpu": task.req_cpus}.get(metric)
                if value:
                    axis.plot(task.rchar / 1024 ** 3, scale(value), "x",
                              color="crimson", markersize=9, markeredgewidth=2)

            axis.set_xscale("log")
            axis.set_yscale(yscale)
            if yscale == "linear":
                axis.set_ylim(bottom=0)
            axis.set_xlabel("task input read (rchar, GB)")
            axis.set_ylabel(ylabel)
            axis.grid(True, which="both", alpha=0.15)

        axes[0].legend(fontsize=7, loc="best", framealpha=0.85)
        fig.suptitle(f"{stats.name}   [{stats.tier}]   "
                     f"{stats.n_tasks} tasks over {stats.n_runs} run(s)   {stats.basis}",
                     fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        path = os.path.join(figures_dir, f"process_{stats.name}.png")
        fig.savefig(path, dpi=110)
        plt.close(fig)
        figures[f"process_{stats.name}"] = path

    return figures


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

HTML_STYLE = """
body { font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
       margin: 0 auto; max-width: 1400px; padding: 2rem 1.5rem; color: #24292f;
       background: #ffffff; line-height: 1.5; }
h1 { font-size: 1.7rem; margin-bottom: .2rem; }
h2 { font-size: 1.25rem; margin-top: 2.2rem; border-bottom: 1px solid #d8dee4;
     padding-bottom: .3rem; }
h3 { font-size: 1rem; margin-top: 1.6rem; color: #444; }
.sub { color: #57606a; font-size: .9rem; margin-top: 0; }
table { border-collapse: collapse; width: 100%; font-size: .82rem; margin-top: .8rem; }
th, td { border: 1px solid #d8dee4; padding: .35rem .5rem; text-align: right; }
th { background: #f6f8fa; text-align: right; position: sticky; top: 0; }
td.name, th.name { text-align: left; font-family: ui-monospace, Menlo, Consolas, monospace; }
tr:nth-child(even) td { background: #fbfcfd; }
.bad { color: #b42318; font-weight: 600; }
.warn { color: #b54708; }
.good { color: #067647; }
img { max-width: 100%; height: auto; display: block; margin: .6rem 0 1.4rem; }
pre { background: #f6f8fa; border: 1px solid #d8dee4; border-radius: 6px;
      padding: .8rem 1rem; overflow-x: auto; font-size: .8rem; }
.note { background: #fff8c5; border: 1px solid #eed888; border-radius: 6px;
        padding: .6rem .9rem; margin: .8rem 0; font-size: .87rem; }
.scroll { overflow-x: auto; }
"""


def _html_escape(text: Any) -> str:
    return (str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def render_html(ordered: Sequence[ProcessStats],
                runs: Sequence[RunContext],
                figures: Dict[str, str],
                config_text: str,
                unmapped: Sequence[str],
                parse_results: Sequence[TraceParseResult],
                cap_findings: Sequence[CapFinding],
                inferred: Sequence[str]) -> str:
    """Build the single self-contained HTML report.

    Written directly rather than through ``bin/dashboard_report.html``: that
    template is a large Plotly application with a per-cell data model that has
    nothing in common with this report, so reusing it would mean maintaining a
    second copy of its placeholder machinery for no gain.
    """
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    n_completed = sum(len(result.completed) for result in parse_results)
    n_failed = sum(len(result.failures) for result in parse_results)
    n_cached = sum(result.n_cached for result in parse_results)

    parts: List[str] = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>Resource efficiency report</title>",
        f"<style>{HTML_STYLE}</style></head><body>",
        "<h1>Cluster resource efficiency</h1>",
        f"<p class='sub'>Generated {stamp} from {len(runs)} run(s) &middot; "
        f"{n_completed} completed task(s), {n_failed} failed, {n_cached} cached rows skipped.</p>",
    ]

    # ── Runs ────────────────────────────────────────────────────────────────
    parts.append("<h2>Runs analysed</h2><div class='scroll'><table>")
    parts.append("<tr><th class='name'>run key</th><th class='name'>trace</th>"
                 "<th>samples</th><th class='name'>protocol</th>"
                 "<th class='name'>mapper</th></tr>")
    for run in runs:
        parts.append(
            f"<tr><td class='name'>{_html_escape(run.key)}</td>"
            f"<td class='name'>{_html_escape(run.trace_path)}</td>"
            f"<td>{len(run.samples)}</td>"
            f"<td class='name'>{_html_escape(run.params.get('protocol', '-'))}</td>"
            f"<td class='name'>{_html_escape(run.params.get('mapping_software', '-'))}</td></tr>")
    parts.append("</table></div>")

    if len(runs) < MIN_FIT_RUNS:
        parts.append(
            "<div class='note'><strong>Single run.</strong> Scatter points still span the "
            "samples in this run, but no scaling exponent is fitted across runs, so every "
            "recommendation is a flat statistic over what was observed. Point the tool at "
            "two or more runs on differently sized datasets to get a scaling fit.</div>")

    missing = sorted({field for result in parse_results for field in result.missing_fields})
    if missing:
        detail = (
            f"<strong>This trace was written with Nextflow's default field set</strong>, "
            f"which records only what a task <em>used</em> &mdash; it is missing "
            f"<code>{_html_escape(', '.join(missing))}</code>. Most likely the run "
            f"predates the <code>trace.fields</code> list now in "
            f"<code>nextflow.config</code>, or was launched with a bare "
            f"<code>-with-trace</code>.")
        if inferred:
            detail += (
                f" The <strong>declared <code>conf/base.config</code> values were "
                f"substituted</strong> for {_html_escape(', '.join(inferred))} so the "
                f"efficiencies below can still be shown. That substitution cannot see "
                f"retry escalation or clamping by the site profile, so treat these "
                f"percentages as close estimates rather than measurements. Re-run the "
                f"pipeline to record the real requests.")
        else:
            detail += " Metrics depending on them are blank."
        parts.append(f"<div class='note'>{detail}</div>")

    if unmapped:
        parts.append(
            f"<div class='note'><strong>Unmapped processes:</strong> "
            f"{_html_escape(', '.join(unmapped))}. No <code>label</code> could be found for "
            f"these in <code>modules/**/main.nf</code>, so they are reported but never "
            f"folded into a tier recommendation.</div>")

    clamped = [stats for stats in ordered if stats.clamped]
    if clamped:
        parts.append(
            f"<div class='note'><strong>Clamped to site limits:</strong> "
            f"{_html_escape(', '.join(stats.name for stats in clamped))}. The recommendation "
            f"exceeded a <code>--max-*</code> ceiling and was capped.</div>")

    # ── Tool-internal memory caps ───────────────────────────────────────────
    notable = [finding for finding in cap_findings if finding.severity != "ok"]
    if notable:
        parts.append("<h2>Tool-internal memory ceilings</h2>")
        parts.append("<p class='sub'>A <code>memory</code> directive only sets what the "
                     "scheduler reserves. These tools take their own ceiling on the command "
                     "line, and it does not move when the directive does &mdash; so the "
                     "recommendations above are only actionable once these agree.</p>")
        parts.append("<div class='scroll'><table><tr><th class='name'>process</th>"
                     "<th class='name'>parameter</th><th>pinned at</th>"
                     "<th class='name'>what to do</th></tr>")
        for finding in notable:
            values = [value for value in finding.pinned_values.values() if value]
            css = {"danger": "bad", "warn": "warn"}.get(finding.severity, "")
            parts.append(
                f"<tr><td class='name'>{_html_escape(finding.process)}</td>"
                f"<td class='name'>{_html_escape(finding.param)}</td>"
                f"<td class='{css}'>{fmt_bytes(max(values)) if values else 'auto'}</td>"
                f"<td class='name' style='text-align:left'>"
                f"{_html_escape(finding.message)}</td></tr>")
        parts.append("</table></div>")

    # ── Overviews ───────────────────────────────────────────────────────────
    for key, title in (("overview_efficiency", "Efficiency by process"),
                       ("overview_tiers", "Spread within each label tier")):
        uri = encode_image(figures.get(key))
        if uri:
            parts.append(f"<h2>{title}</h2><img src='{uri}' alt='{title}'>")

    # ── Recommendation table ────────────────────────────────────────────────
    parts.append("<h2>Per-process summary and recommendations</h2><div class='scroll'><table>")
    parts.append(
        "<tr><th class='name'>process</th><th class='name'>label</th><th>runs</th>"
        "<th>tasks</th><th>retries</th><th>max peak RSS</th><th>requested</th>"
        "<th>mem eff</th><th>cpu eff</th><th>p95 runtime</th><th>time eff</th>"
        "<th>b</th><th>rec cpus</th><th>rec memory</th><th>rec time</th>"
        "<th class='name'>basis</th></tr>")
    for stats in ordered:
        fit = stats.fits.get("memory")
        exponent = f"{fit.slope:.2f}" if fit and fit.trusted else "&ndash;"
        efficiency = stats.memory_efficiency
        css = "bad" if (efficiency and efficiency > 0.9) else (
            "warn" if (efficiency is not None and efficiency < 0.25) else "good")
        parts.append(
            f"<tr><td class='name'>{_html_escape(stats.name)}</td>"
            f"<td class='name'>{_html_escape(stats.tier)}</td>"
            f"<td>{stats.n_runs}</td><td>{stats.n_tasks}</td>"
            f"<td>{stats.n_retries or ''}</td>"
            f"<td>{fmt_bytes(stats.max_peak_rss)}</td>"
            f"<td>{fmt_bytes(stats.current_memory)}</td>"
            f"<td class='{css}'>{fmt_pct(efficiency)}</td>"
            f"<td>{fmt_pct(stats.cpu_efficiency)}</td>"
            f"<td>{fmt_duration(stats.p95_realtime)}</td>"
            f"<td>{fmt_pct(stats.time_efficiency)}</td>"
            f"<td>{exponent}</td>"
            f"<td>{stats.rec_cpus or ''}</td>"
            f"<td>{fmt_bytes(stats.rec_memory) if stats.rec_memory else ''}</td>"
            f"<td>{fmt_duration(stats.rec_time) if stats.rec_time else ''}</td>"
            f"<td class='name'>{_html_escape(stats.basis)}</td></tr>")
    parts.append("</table></div>")

    # ── Failures ────────────────────────────────────────────────────────────
    killed = [(stats, task) for stats in ordered for task in stats.failures if task.killed]
    if killed:
        parts.append("<h2>Tasks the scheduler killed</h2>")
        parts.append("<p class='sub'>These are the strongest signal in the report: the "
                     "request shown was demonstrably too small. Recommendations for these "
                     "processes were raised regardless of what the surviving tasks used.</p>")
        parts.append("<div class='scroll'><table><tr><th class='name'>process</th>"
                     "<th class='name'>run</th><th class='name'>tag</th><th>exit</th>"
                     "<th>attempt</th><th>requested memory</th><th>requested time</th>"
                     "<th class='name'>job id</th></tr>")
        for stats, task in killed:
            parts.append(
                f"<tr><td class='name'>{_html_escape(stats.name)}</td>"
                f"<td class='name'>{_html_escape(task.run_key)}</td>"
                f"<td class='name'>{_html_escape(task.tag)}</td>"
                f"<td class='bad'>{task.exit_code}</td><td>{task.attempt}</td>"
                f"<td>{fmt_bytes(task.req_memory)}</td>"
                f"<td>{fmt_duration(task.req_time)}</td>"
                f"<td class='name'>{_html_escape(task.native_id)}</td></tr>")
        parts.append("</table></div>")

    # ── Per-process panels ──────────────────────────────────────────────────
    panel_keys = [key for key in figures if key.startswith("process_")]
    if panel_keys:
        parts.append("<h2>Usage against dataset size, per process</h2>")
        parts.append("<p class='sub'>One point per task; x is the bytes that task read "
                     "(<code>rchar</code>). Colours are runs. The grey line is the current "
                     "request, green the recommendation, red crosses are killed tasks. A "
                     "fitted line is drawn only where the fit passed every quality gate.</p>")
        for stats in ordered:
            uri = encode_image(figures.get(f"process_{stats.name}"))
            if uri:
                parts.append(f"<h3>{_html_escape(stats.name)}</h3>")
                parts.append(f"<img src='{uri}' alt='{_html_escape(stats.name)}'>")

    # ── Emitted config ──────────────────────────────────────────────────────
    parts.append("<h2>Proposed Nextflow configuration</h2>")
    parts.append("<p class='sub'>Save as <code>conf/resources_tuned.config</code> (or re-run "
                 "with <code>--apply</code>) and it will be picked up automatically by "
                 "<code>submit_nextflow.sh</code>.</p>")
    parts.append(f"<pre>{_html_escape(config_text)}</pre>")

    parts.append("</body></html>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Report HPC resource efficiency across Nextflow runs and "
                    "propose retuned process resources.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # ── Inputs ──────────────────────────────────────────────────────────────
    parser.add_argument("-r", "--results", action="append", default=[], metavar="DIR",
                        help="Pipeline result directory (a params.outdir) holding "
                             "pipeline_info/. Repeatable; more runs give a better "
                             "scaling fit.")
    parser.add_argument("-R", "--results-root", metavar="DIR",
                        help="Directory whose immediate subdirectories are result "
                             "directories. Convenient when every run went to its own outdir.")
    parser.add_argument("-p", "--pipeline-dir", metavar="DIR",
                        help="Pipeline checkout holding modules/ and conf/base.config, used "
                             "to recover the process->label mapping. Default: the parent of "
                             "this script.")

    # ── Outputs ─────────────────────────────────────────────────────────────
    parser.add_argument("-o", "--output", metavar="DIR",
                        help="Report directory. Default: ./resource_efficiency_<timestamp>/")
    parser.add_argument("--json", metavar="PATH",
                        help="Also write the full analysis as JSON.")
    parser.add_argument("--emit-config", metavar="PATH",
                        help="Write the tuned Nextflow config to this path.")
    parser.add_argument("--apply", action="store_true",
                        help="Write the tuned config to <pipeline-dir>/conf/resources_tuned.config, "
                             "where submit_nextflow.sh picks it up on the next run.")
    parser.add_argument("--dynamic", action="store_true",
                        help="Also write input-size dependent `memory` directives for the "
                             "processes whose usage provably tracks their input, as "
                             "paste-ready snippets for the module files. A flat request "
                             "has to cover the largest dataset; a closure over the input "
                             "size does not.")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip the figures and the HTML report. Needs no matplotlib.")
    parser.add_argument("--max-panels", type=int, default=40, metavar="N",
                        help="Cap on per-process figures, most expensive first. Default: 40.")

    # ── Tuning policy ───────────────────────────────────────────────────────
    parser.add_argument("--target-size", type=float, metavar="GB",
                        help="Extrapolate recommendations to a dataset of this size, using "
                             "the fitted scaling exponent where the fit is trustworthy.")
    parser.add_argument("--min-tasks", type=int, default=3, metavar="N",
                        help="Observations required before a process is retuned. Default: 3.")
    parser.add_argument("--safety-memory", type=float, default=DEFAULT_SAFETY_MEMORY,
                        metavar="F",
                        help=f"Factor applied to the observed memory maximum. "
                             f"Default: {DEFAULT_SAFETY_MEMORY}.")
    parser.add_argument("--safety-time", type=float, default=DEFAULT_SAFETY_TIME, metavar="F",
                        help=f"Factor applied to p95 runtime. Default: {DEFAULT_SAFETY_TIME}.")

    # ── Site limits ─────────────────────────────────────────────────────────
    parser.add_argument("--max-cpus", type=int, metavar="N",
                        help="Cap recommended cpus. The CRG limits live in the external "
                             "nf-core/configs crg profile and will clamp anything larger anyway.")
    parser.add_argument("--max-memory", type=float, metavar="GB",
                        help="Cap recommended memory, in GB.")
    parser.add_argument("--max-time", type=float, metavar="H",
                        help="Cap recommended walltime, in hours.")

    args = parser.parse_args(argv)

    if not args.results and not args.results_root:
        parser.error("give at least one --results DIR (or --results-root DIR)")

    if not args.pipeline_dir:
        args.pipeline_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    if not args.output:
        stamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        args.output = os.path.join(os.getcwd(), f"resource_efficiency_{stamp}")

    return args


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    # ── Discover runs ───────────────────────────────────────────────────────
    runs = discover_runs(args.results, args.results_root)
    if not runs:
        sys.stderr.write(
            "Error: no execution traces found.\n"
            "  Each run writes <outdir>/pipeline_info/execution_trace_<timestamp>.txt.\n"
            "  If that file is missing, the run predates the trace block in nextflow.config;\n"
            "  re-run the pipeline and try again.\n")
        sys.exit(1)

    print(f"Found {len(runs)} run(s):")
    for run in runs:
        print(f"  {run.describe()}")

    # ── Label mapping ───────────────────────────────────────────────────────
    module_labels = scan_module_labels(args.pipeline_dir)
    tiers = parse_base_config(os.path.join(args.pipeline_dir, "conf", "base.config"))

    # Discover the tier names from base.config before any ProcessLabels.tier is read,
    # so a tier added there is picked up without touching this file.
    adopt_resource_labels(tiers)

    if not module_labels:
        sys.stderr.write(f"Warning: no processes found under {args.pipeline_dir}/modules; "
                         f"label mapping will be empty. Check --pipeline-dir.\n")

    # A label no withLabel block defines means this run is measured against the wrong
    # tier, silently. Say so rather than folding the process into __default__.
    unrecognised = sorted({label
                           for info in module_labels.values()
                           for label in info.unrecognised})
    if unrecognised:
        sys.stderr.write(
            f"Warning: {', '.join(unrecognised)} used by a module but not declared in "
            f"conf/base.config; processes carrying only those labels are treated as "
            f"unlabelled.\n")

    # ── Parse traces ────────────────────────────────────────────────────────
    parse_results: List[TraceParseResult] = []
    all_tasks: List[TaskRecord] = []
    all_failures: List[TaskRecord] = []
    unmapped: set = set()

    for run in runs:
        result = parse_trace(run.trace_path, run.key)
        parse_results.append(result)
        for task in result.completed + result.failures:
            name = simple_process_name(task.process, task.name)
            module_process, tier = resolve_label(name, module_labels)
            if tier is None:
                unmapped.add(name)
            task.module_process = module_process
            task.tier = tier or DEFAULT_TIER
        all_tasks.extend(result.completed)
        all_failures.extend(result.failures)

    if not all_tasks:
        sys.stderr.write("Error: the traces contain no COMPLETED tasks to analyse.\n")
        sys.exit(1)

    print(f"\nParsed {len(all_tasks)} completed task(s), {len(all_failures)} failed, "
          f"{sum(r.n_cached for r in parse_results)} cached row(s) skipped.")
    if unmapped:
        sys.stderr.write(f"Warning: no label found for {len(unmapped)} process(es): "
                         f"{', '.join(sorted(unmapped))}\n")

    # ── Backfill requested resources the trace did not record ───────────────
    # Without a denominator every efficiency figure, and both overview charts,
    # would be blank. The declared tier value is the best stand-in available; see
    # backfill_requests() for where it is only an approximation.
    inferred = backfill_requests(list(all_tasks) + list(all_failures), tiers)
    if inferred:
        print(f"\nNote: this trace does not record {', '.join(inferred)}. "
              f"Those are Nextflow's\n"
              f"      defaults being used rather than the field list in "
              f"nextflow.config -- most likely\n"
              f"      the run predates it, or used a bare `-with-trace`. The declared "
              f"conf/base.config\n"
              f"      values were substituted so efficiencies can still be shown; they "
              f"do not account\n"
              f"      for retry escalation or clamping by the site profile. Re-run to "
              f"record them properly.")

    # ── Aggregate and recommend ─────────────────────────────────────────────
    stats_by_name = aggregate(all_tasks, all_failures)
    for stats in stats_by_name.values():
        recommend(stats, args)

    ordered = sorted(stats_by_name.values(), key=lambda s: -s.cpu_hours)
    cap_findings = check_tool_caps(stats_by_name, runs)

    # ── Emit ────────────────────────────────────────────────────────────────
    os.makedirs(args.output, exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    print()
    print(render_table(ordered))

    notable = [finding for finding in cap_findings if finding.severity != "ok"]
    if notable:
        print()
        print("Tool-internal memory ceilings (a bigger request alone will not move these):")
        for finding in notable:
            marker = {"danger": "!!", "warn": " !", "info": "  "}[finding.severity]
            print(f"  {marker} {finding.process} / {finding.param}")
            print(f"       {finding.message}")

    tsv_path = os.path.join(args.output, f"resource_efficiency_{stamp}.tsv")
    write_tsv(tsv_path, ordered)

    plans = build_tier_plans(stats_by_name, tiers, module_labels, args)
    config_text = render_config(plans, tiers, runs, len(all_tasks), cap_findings, args)
    config_path = os.path.join(args.output, f"resources_tuned_{stamp}.config")
    with open(config_path, "w", encoding="utf-8") as handle:
        handle.write(config_text)

    dynamic_path: Optional[str] = None
    if args.dynamic:
        dynamic_path = os.path.join(args.output, f"resources_dynamic_{stamp}.config")
        with open(dynamic_path, "w", encoding="utf-8") as handle:
            handle.write(render_dynamic(ordered, tiers, args))

    figures: Dict[str, str] = {}
    html_path: Optional[str] = None
    if not args.no_plots:
        figures = make_plots(ordered, runs, os.path.join(args.output, "figures"), args)
        html_path = os.path.join(args.output, f"resource_efficiency_{stamp}.html")
        with open(html_path, "w", encoding="utf-8") as handle:
            handle.write(render_html(ordered, runs, figures, config_text,
                                     sorted(unmapped), parse_results, cap_findings,
                                     inferred))

    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump(build_json(ordered, runs, sorted(unmapped), parse_results, cap_findings),
                      handle, indent=2, default=str)

    if args.emit_config:
        with open(args.emit_config, "w", encoding="utf-8") as handle:
            handle.write(config_text)

    applied_path: Optional[str] = None
    if args.apply:
        applied_path = os.path.join(args.pipeline_dir, "conf", "resources_tuned.config")
        os.makedirs(os.path.dirname(applied_path), exist_ok=True)
        with open(applied_path, "w", encoding="utf-8") as handle:
            handle.write(config_text)

    # ── Summary ─────────────────────────────────────────────────────────────
    print()
    print(f"Report directory : {args.output}")
    print(f"  per-process TSV: {tsv_path}")
    print(f"  tuned config   : {config_path}")
    if dynamic_path:
        print(f"  dynamic sizing : {dynamic_path}")
    if html_path:
        print(f"  HTML report    : {html_path}")
        print(f"  figures        : {len(figures)} PNG(s) under {args.output}/figures")
    if args.json:
        print(f"  JSON           : {args.json}")

    if applied_path:
        print()
        print(f"Applied to {applied_path}.")
        print("submit_nextflow.sh adds it automatically on the next run; to use it directly:")
        print(f"    nextflow run ... -c {applied_path}")
        print("Delete that file to revert to conf/base.config.")
    else:
        print()
        print("Nothing outside the report directory was modified. To use the recommendations:")
        print(f"    cp {config_path} {os.path.join(args.pipeline_dir, 'conf', 'resources_tuned.config')}")
        print("  or re-run with --apply.")

    if len(runs) < MIN_FIT_RUNS:
        print()
        print(f"Note: only {len(runs)} run analysed, so no scaling exponent was fitted and every")
        print("      recommendation is a flat statistic. Point --results at runs on differently")
        print("      sized datasets to predict how resources scale.")


if __name__ == "__main__":
    main()
