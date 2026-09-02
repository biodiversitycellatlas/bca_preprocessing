#!/usr/bin/env python3
"""
Embed data into HTML dashboard summarising a scRNA-seq preprocessing run.

Usage modes
-----------
Pipeline mode (Nextflow):
    generate_dashboard.py --template bin/dashboard_report.html --output /path/to/output.html --samplesheet /path/to/samplesheet.csv
        --analytical_samples ... [--star_logs ...] [--af_meta_info ...] ...

Standalone mode:
    generate_dashboard.py --template bin/dashboard_report.html --output /path/to/output.html --result-dir /path/to/results
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
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Ordered longest-first so "_geneext_starsolo" is checked before "_starsolo".
_ANALYTICAL_SUFFIXES: List[str] = [
    "_geneext_starsolo",
    "_starsolo",
    "_alevinfry",
]

# STARsolo gene-model directories, in preference order.
_SOLO_GENE_DIRS: List[str] = ["GeneFull_Ex50pAS", "Gene"]

# Subdirectory of a gene-model directory holding the second-derivative outputs.
_SECONDDERIV_DIR: str = "filtered_secondderiv"

# Candidate locations for optional cell-filtering outputs (checked in order).
_CELLSWEEP_ROOTS:  List[str] = ["cell_filtering/cellsweep", "cellsweep"]

# Candidate locations for Kraken sankey HTML files.
_SANKEY_ROOTS: List[str] = ["kraken", "taxonomy"]

# Suffixes alevin-fry appends to the spliced / unspliced / ambiguous column blocks of a USA-mode count matrix.
_USA_SUFFIX_RE = re.compile(r"-[SUA]$")

# ── Payload budgets ──────────────────────────────────────────────────────────
# The report embeds every curve it draws as JSON inside the HTML, so its size is
# driven by the number of *barcodes*, not the number of samples: a deeply
# sequenced library has millions of barcodes with at least one UMI, and writing
# one JSON number per barcode per series runs into hundreds of MB.  The budgets
# below cap each series at the point where adding more data stops changing the
# picture -- the barcode-rank and cumulative curves are drawn on log axes, where
# a curve thinned uniformly in log10(rank) is indistinguishable from the full
# one, and the per-cell scatters saturate their pixels long before the full
# ambient cloud is drawn.  Every number that the dashboard *reports* (cell
# counts, medians, thresholds) is computed upstream on the complete data and is
# unaffected by this thinning.
_KNEE_MAX_POINTS: int = 3000        # barcode-rank curve, drawn on log-log axes
_KNEE_HEAD_RANKS: int = 500         # leading ranks kept untouched (the called cells)
_SD_MAX_POINTS:   int = 3000        # second-derivative cumulative-UMI curve
_PERCELL_MAX_NONCELLS: int = 50000  # grey background cloud of the per-cell scatters

# Minimum UMIs per barcode considered by the second-derivative cell calling;
# must match MIN_UMIS in bin/secondderiv_cellcalling.py.
_SD_MIN_UMIS: int = 100

# Per-cell metric columns the report actually plots.  ``Cell`` (the barcode
# string) is dropped from the embedded copy -- nothing reads it, and it is the
# single widest column.  The published *_metrics.csv / *_metrics.json keep it.
_PERCELL_COLUMNS: List[str] = [
    "IntronicPercent", "MTPercent", "rRNAPercent", "TotalReads", "IsCell",
]
_PERCELL_PCT_COLUMNS: List[str] = ["IntronicPercent", "MTPercent", "rRNAPercent"]

# ── GeneExt ──────────────────────────────────────────────────────────────────
# GeneExt writes a standalone HTML report next to the extended GTF with every
# number it computed embedded as a single JSON object.  The dashboard re-uses
# that object rather than re-deriving anything from the GTF, but keeps only the
# run-level summary and the two histograms it redraws: the report's ``ext_table``
# carries one row per extended gene and ``orphan_bed`` the full orphan-peak BED,
# which together dwarf every other block and are already one click away in the
# GeneExt report itself.
_GENEEXT_SUMMARY_KEYS: List[str] = [
    "n_genes", "n_extended", "pct_extended",
    "min_ext", "median_ext", "mean_ext", "max_ext",
    "n_genic_peaks", "n_noov_peaks", "n_orphan_peaks", "orphan_warning",
    "cov_percentile", "cov_threshold", "n_reads", "subsampled",
    "input_file", "output_file", "genome_fixed", "run_date", "run_args",
]

# Blocks copied verbatim from the GeneExt payload (all small, fixed-size).
_GENEEXT_BLOCK_KEYS: List[str] = ["ext_hist", "cov_hist", "peak_flow", "log_notes"]


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse CLI arguments.

    Both ``--analytical_samples`` / ``--samplesheet`` (pipeline mode) and
    ``--result-dir`` (standalone mode) are optional here; mutual-exclusion
    validation is performed in :func:`main`.
    """
    parser = argparse.ArgumentParser(
        description="Generate HTML Dashboard for scRNA-seq preprocessing results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # ── Required for both modes ──────────────────────────────────────────────
    parser.add_argument("--template", required=True, help="Path to HTML template")
    parser.add_argument("--output",   required=True, help="Output HTML filename")

    # ── Core metadata (optional overrides; auto-detected in standalone mode) ─
    parser.add_argument("--project",  default="Biodiversity Cell Atlas", help="Project name")
    parser.add_argument("--pipeline", default="bca-preprocessing",       help="Pipeline name")
    parser.add_argument("--version",  default=None, help="Pipeline version (auto-detected from params JSON when omitted)")
    parser.add_argument("--commit",   default=None, help="Git commit hash  (auto-detected from params JSON when omitted)")

    # ── Standalone mode ──────────────────────────────────────────────────────
    parser.add_argument(
        "--result-dir", dest="result_dir", default=None,
        help="Root of the pipeline result directory (standalone mode). "
             "When supplied, all file-list arguments are ignored.",
    )

    # ── Pipeline mode: configurations & manifest ─────────────────────────────
    parser.add_argument("--run_config",         default=None, help="Path to run_config_date.txt")
    parser.add_argument("--samplesheet",        default=None, help="Path to samplesheet.csv")
    parser.add_argument("--analytical_samples", default=None, help="CSV manifest from Nextflow")

    # ── Pipeline mode: STARsolo file lists ───────────────────────────────────
    parser.add_argument("--star_logs",      nargs="*", help="STAR Log.final.out files")
    parser.add_argument("--star_summaries", nargs="*", help="STARsolo Summary.csv files")
    parser.add_argument("--star_full_logs", nargs="*", help="STAR Log.out files (UMI threshold)")
    parser.add_argument("--cell_stats",     nargs="*", help="STARsolo CellReads.stats files")

    # ── Pipeline mode: alevin-fry file lists ─────────────────────────────────
    parser.add_argument("--af_meta_info",  nargs="*", help="meta_info.json files")
    parser.add_argument("--af_quant_json", nargs="*", help="quant.json files")
    parser.add_argument("--af_cell_meta",  nargs="*", help="cell_meta.tsv files")
    parser.add_argument("--af_mat_cols",   nargs="*", help="quants_mat_cols.txt files (authoritative gene count)")

    # ── Pipeline mode: visualisation file lists ──────────────────────────────
    parser.add_argument("--sankey_files",    nargs="*", help="*_kraken.sankey.html files")
    parser.add_argument("--per_cell_files",  nargs="*", help="Per-cell JSON files")
    parser.add_argument("--saturation_imgs", nargs="*", help="Saturation PNG images")
    parser.add_argument("--residuals_imgs",  nargs="*", help="Residuals PNG images")
    parser.add_argument("--saturation_logs",nargs="*", help="saturation.log files")
    parser.add_argument("--mt_rrna_metrics",nargs="*", help="*_mt_rrna_metrics.txt files")
    parser.add_argument("--knee_files",      nargs="*", help="UMIperCellSorted.txt files")
    parser.add_argument("--secondderiv_knee",  nargs="*", help="*_knee_data.json files from the second-derivative cell calling")
    parser.add_argument("--secondderiv_stats", nargs="*", help="*_secondderiv_statistics.json files from the second-derivative matrix filtering")

    parser.add_argument("--cellsweep_tables",        nargs="*")
    parser.add_argument("--cellsweep_plots_contrib", nargs="*")
    parser.add_argument("--cellsweep_plots_umap",    nargs="*")

    # ── Pipeline mode: GeneExt (run-level, not per sample) ───────────────────
    parser.add_argument("--geneext_report", nargs="*", help="GeneExt *.Report.html file")
    parser.add_argument("--geneext_log",    nargs="*", help="GeneExt *.GeneExt.log file")

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Formatters & encoders
# ---------------------------------------------------------------------------

def to_pct(val: Any, precision: int = 2) -> str:
    """Convert a fraction (0–1) to a percentage string, or return ``"N/A"``."""
    if val is None or str(val).strip().upper() == "N/A" or val == "":
        return "N/A"
    try:
        return f"{float(val) * 100:.{precision}f}%"
    except (ValueError, TypeError):
        return "N/A"


def encode_image(image_path: Optional[str]) -> Optional[str]:
    """Return a PNG data-URI for HTML embedding, or ``None`` if unreadable."""
    if not image_path or not os.path.exists(image_path):
        return None
    try:
        with open(image_path, "rb") as fh:
            encoded = base64.b64encode(fh.read()).decode("utf-8")
        return f"data:image/png;base64,{encoded}"
    except Exception as exc:
        sys.stderr.write(f"Warning: could not encode image {image_path}: {exc}\n")
        return None


def safe_float(val: Any) -> Optional[float]:
    """Return ``float(val)`` or ``None`` on failure."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def fmt(val: Any) -> str:
    """Format integers with comma separators; pass other values through as-is."""
    if val is None or val == "":
        return ""
    try:
        s = str(val).strip()
        if s.isdigit():
            return f"{int(s):,}"
    except Exception:
        pass
    return str(val)


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def extract_sankey_data(html_path: Optional[str]) -> Optional[dict]:
    """Extract the htmlwidgets JSON payload from a sankey HTML file.

    The widget ID is whatever htmlwidgets generated -- older versions emit a hex
    digest, newer ones an arbitrary alphanumeric string -- and the attribute
    order varies between renderers, so the payload is located by its
    ``data-for="htmlwidget-…"`` marker rather than by an exact tag spelling.
    """
    if not html_path or not os.path.exists(html_path):
        return None
    try:
        with open(html_path, "r") as fh:
            content = fh.read()
        match = re.search(
            r'<script[^>]*\bdata-for="htmlwidget-[^"]+"[^>]*>(.*?)</script>',
            content, re.DOTALL,
        )
        if not match:
            # Attribute order is not fixed: type may follow data-for instead.
            match = re.search(
                r'<script[^>]*\btype="application/json"[^>]*>(.*?)</script>',
                content, re.DOTALL,
            )
        if match:
            return json.loads(match.group(1)).get("x")
        sys.stderr.write(
            f"Warning: no htmlwidgets payload found in {html_path}; "
            "the taxonomy tab will be empty for this sample\n"
        )
    except Exception as exc:
        sys.stderr.write(f"Warning: could not read sankey {html_path}: {exc}\n")
    return None


def parse_saturation_log(path: Optional[str]) -> str:
    """Return reads needed for target saturation (e.g. ``'22.5 M'``), or ``'N/A'``."""
    if not path or not os.path.exists(path):
        return "N/A"
    try:
        with open(path, "r") as fh:
            content = fh.read()
        match = re.search(
            r"To achieve a saturation of [0-9.]+, ninput should be approximately:\s*(.*?)\s*reads",
            content,
        )
        if match:
            return match.group(1).strip()
    except Exception:
        pass
    return "N/A"


def extract_geneext_report_data(html_path: Optional[str]) -> Dict[str, Any]:
    """Extract the embedded ``const D = {...}`` payload from a GeneExt report.

    The object is written on a single line by GeneExt, but it contains braces and
    escaped quotes throughout, so it is decoded with :class:`json.JSONDecoder`
    from the first ``{`` after the declaration rather than matched with a regex.

    Returns ``{}`` when the file is missing or holds no recognisable payload.
    """
    if not html_path or not os.path.exists(html_path):
        return {}
    try:
        with open(html_path, "r", encoding="utf-8", errors="replace") as fh:
            content = fh.read()
        match = re.search(r"\bconst\s+D\s*=\s*", content)
        if not match:
            sys.stderr.write(
                f"Warning: no GeneExt data payload found in {html_path}; "
                "falling back to the GeneExt log\n"
            )
            return {}
        obj, _ = json.JSONDecoder().raw_decode(content, match.end())
        return obj if isinstance(obj, dict) else {}
    except Exception as exc:
        sys.stderr.write(f"Warning: could not read GeneExt report {html_path}: {exc}\n")
        return {}


def parse_geneext_log(path: Optional[str]) -> Dict[str, Any]:
    """Recover GeneExt's headline numbers from its plain-text log.

    Used when the HTML report is absent -- a GeneExt older than the one that
    writes it -- so the dashboard can still report what the extension did.  Only
    the three lines GeneExt prints as its result are read; anything the report
    would have added (the histograms, the peak-filtering flow) has no equivalent
    in the log and is simply left out.
    """
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            content = fh.read()
    except Exception as exc:
        sys.stderr.write(f"Warning: could not read GeneExt log {path}: {exc}\n")
        return {}

    summary: Dict[str, Any] = {}
    peak_flow: Dict[str, Any] = {}

    extended = re.search(r"Extended\s+(\d+)\s*/\s*(\d+)\s+genes", content)
    if extended:
        n_extended, n_genes = int(extended.group(1)), int(extended.group(2))
        summary["n_extended"] = n_extended
        summary["n_genes"]    = n_genes
        if n_genes:
            summary["pct_extended"] = round(100.0 * n_extended / n_genes, 1)

    median = re.search(r"Median extension length:\s*([\d.]+)\s*bp", content)
    if median:
        summary["median_ext"] = float(median.group(1))

    retained = re.search(
        r"Retained\s+(\d+)\s*/\s*(\d+)\s+\(\s*[\d.]+\s*%\s*\)\s+intergenic peaks", content
    )
    if retained:
        peak_flow["passed_filtering"] = int(retained.group(1))
        summary["n_noov_peaks"]       = int(retained.group(2))
    if "n_extended" in summary:
        peak_flow["assigned_to_genes"] = summary["n_extended"]

    notes = [
        line.strip() for line in content.splitlines()
        if "warning" in line.lower() and line.strip()
    ]

    if not summary and not peak_flow:
        return {}
    return {"summary": summary, "peak_flow": peak_flow, "log_notes": notes}


def build_geneext_payload(
    report_path: Optional[str],
    log_path: Optional[str],
) -> Dict[str, Any]:
    """Assemble the gene-extension block the dashboard embeds.

    Prefers GeneExt's own HTML report, which carries every statistic it computed,
    and falls back to the log when the report is unavailable.  ``source`` records
    which of the two the numbers came from, so the tab can say when it is showing
    the reduced, log-derived view.

    Returns ``{}`` when neither file yields anything, which is what hides the tab.
    """
    payload: Dict[str, Any] = {}

    data = extract_geneext_report_data(report_path)
    if data:
        summary = data.get("summary") or {}
        payload["summary"] = {
            k: summary[k] for k in _GENEEXT_SUMMARY_KEYS if k in summary
        }
        for key in _GENEEXT_BLOCK_KEYS:
            if data.get(key):
                payload[key] = data[key]

        fix_info = data.get("fix_info") or {}
        if fix_info.get("extension_param"):
            payload["extension_param"] = fix_info["extension_param"]
        if fix_info.get("steps"):
            # Only the steps that actually ran are worth a line in the report
            payload["genome_fix"] = {
                name: step for name, step in fix_info["steps"].items()
                if isinstance(step, dict) and step.get("applied")
            }
        if data.get("mapping_stats"):
            payload["mapping_stats"] = data["mapping_stats"]
        payload["source"] = "report"

    else:
        from_log = parse_geneext_log(log_path)
        if not from_log:
            return {}
        payload = {k: v for k, v in from_log.items() if v}
        payload["source"] = "log"

    if report_path:
        payload["report_file"] = os.path.basename(report_path)
    if log_path:
        payload["log_file"] = os.path.basename(log_path)
    return payload


def parse_run_config(path: Optional[str]) -> Dict[str, str]:
    """Parse a ``KEY=VALUE`` text file into a dict; return ``{}`` on failure."""
    config: Dict[str, str] = {}
    if not path or not os.path.exists(path):
        return config
    try:
        with open(path, "r") as fh:
            for line in fh:
                if "=" in line:
                    key, val = line.split("=", 1)
                    config[key.strip()] = val.strip()
    except Exception:
        pass
    return config


def parse_cell_filtering_table(path: Optional[str]) -> List[Dict[str, Any]]:
    """Parse a cellsweep output table into a list of row dicts."""
    if not path or not os.path.exists(path):
        return []
    data: List[Dict[str, Any]] = []
    try:
        with open(path, "r", newline="", encoding="utf-8-sig") as fh:
            header_line = fh.readline()
            delimiter = "\t" if "\t" in header_line else ","
            fh.seek(0)
            for row in csv.DictReader(fh, delimiter=delimiter):
                cleaned = {k.strip(): v.strip() for k, v in row.items() if k}
                data.append({
                    "gene":        cleaned.get("gene", "N/A"),
                    "ambient_hat": safe_float(cleaned.get("ambient_hat")),
                })
        data.sort(key=lambda x: x["ambient_hat"] if x["ambient_hat"] is not None else -1, reverse=True)
        return data[:100]
    except Exception:
        return []


def parse_star_log(path: Optional[str]) -> Dict[str, str]:
    """Parse a STAR ``Log.final.out`` (``key | value``) into a dict."""
    data: Dict[str, str] = {}
    if not path or not os.path.exists(path):
        return data
    try:
        with open(path, "r") as fh:
            for line in fh:
                if "|" in line:
                    parts = line.split("|", 1)
                    data[parts[0].strip()] = parts[1].strip()
    except Exception:
        pass
    return data


def parse_starsolo_summary(path: Optional[str]) -> Dict[str, str]:
    """Parse a STARsolo ``Summary.csv`` into a dict."""
    data: Dict[str, str] = {}
    if not path or not os.path.exists(path):
        return data
    try:
        with open(path, "r") as fh:
            for row in csv.reader(fh):
                if len(row) >= 2:
                    data[row[0].strip()] = row[1].strip()
    except Exception:
        pass
    return data


def parse_mt_rrna_metrics(path: Optional[str]) -> Dict[str, str]:
    """Parse an MT/rRNA metrics file (``key,value`` lines) into a dict."""
    data: Dict[str, str] = {}
    if not path or not os.path.exists(path):
        return data
    try:
        with open(path, "r") as fh:
            for line in fh:
                if "," in line:
                    key, val = line.strip().split(",", 1)
                    data[key.strip()] = val.strip()
    except Exception:
        pass
    return data


def parse_starsolo_intronic(path: Optional[str]) -> Any:
    """Compute intronic fraction from a STARsolo ``CellReads.stats`` file.

    STARsolo records the region flags (``exonic``, ``intronic``, ``exonicAS``,
    ``intronicAS``, ``mito``) only for reads with a unique genomic alignment, so
    ``genomeU`` is the denominator; ``genomeM`` reads can never appear in the
    numerator and would deflate the fraction.  The ``CBnotInPasslist`` summary
    row aggregates reads whose barcode was not in the passlist and is skipped.

    Columns are resolved by header name rather than position, since the column
    set varies between STAR versions (some emit ``cbPerfectU`` / ``cbMMuniqueU``
    / ``cbMMmultipleU`` after ``featureM``).

    The fraction is computed over all barcodes, not only the called cells, so it
    is a library-level rate.

    Returns the fraction as a ``float``, or ``"N/A"`` when unavailable.
    """
    if not path or not os.path.exists(path):
        return "N/A"
    try:
        with open(path, "r") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {name: i for i, name in enumerate(header)}
            if "intronic" not in idx or "genomeU" not in idx:
                sys.stderr.write(
                    f"Warning: 'intronic'/'genomeU' columns not found in {path}\n"
                )
                return "N/A"

            i_intronic, i_genome_u = idx["intronic"], idx["genomeU"]
            needed = max(i_intronic, i_genome_u)

            intronic_sum = 0
            genome_unique = 0
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= needed or parts[0] == "CBnotInPasslist":
                    continue
                intronic_sum += int(parts[i_intronic])
                genome_unique += int(parts[i_genome_u])

        if genome_unique > 0:
            return intronic_sum / genome_unique
    except Exception:
        pass
    return "N/A"


def extract_umi_cutoff(log_out_path: Optional[str], cfg_name: str = "GeneFull_Ex50pAS") -> int:
    """Extract ``nUMImin`` from a STAR ``Log.out`` Solo post-map block.

    Returns ``0`` when the value cannot be determined.
    """
    if not log_out_path or not os.path.exists(log_out_path):
        return 0
    try:
        in_block = False
        with open(log_out_path, "r") as fh:
            for line in fh:
                if f"Starting Solo post-map for {cfg_name}" in line:
                    in_block = True
                elif in_block and "Starting Solo post-map for" in line:
                    break
                elif in_block and "cellFiltering" in line:
                    match = re.search(r"nUMImin=(\d+)", line)
                    if match:
                        return int(match.group(1))
    except Exception:
        pass
    return 0


def count_alevin_genes(quant: Dict[str, Any], cols_path: Optional[str]) -> Any:
    """Number of distinct genes in an alevin-fry reference.

    ``quant.json``'s ``num_genes`` counts *matrix columns*.  In USA mode -- which
    this pipeline always runs, since ``alevin-fry quant`` is given a 3-column
    ``t2g_3col.tsv`` -- the matrix carries three blocks per gene (spliced,
    unspliced, ambiguous), so that value is three times the gene count.

    ``quants_mat_cols.txt`` names those columns and is authoritative, so it is
    preferred when available.  The USA collapse is applied only when stripping
    the ``-S``/``-U``/``-A`` suffixes yields exactly one third as many distinct
    names, so a non-USA reference -- or a gene legitimately ending in ``-S`` --
    is never mis-collapsed.  ``quant.json`` is the fallback.

    Note this is the size of the reference, not the number of genes with
    non-zero counts; the dashboard labels it accordingly.
    """
    if cols_path and os.path.exists(cols_path):
        try:
            with open(cols_path, "r") as fh:
                names = [line.strip() for line in fh if line.strip()]
            if names:
                if len(names) % 3 == 0:
                    stripped = {_USA_SUFFIX_RE.sub("", n) for n in names}
                    if len(stripped) == len(names) // 3:
                        return len(stripped)
                return len(set(names))
        except Exception:
            sys.stderr.write(f"Warning: could not read {cols_path}\n")

    n_cols = quant.get("num_genes", quant.get("num_targets"))
    try:
        n_cols = int(n_cols)
    except (TypeError, ValueError):
        return "N/A"
    if quant.get("usa_mode") and n_cols % 3 == 0:
        return n_cols // 3
    return n_cols


def manifest_version() -> Optional[str]:
    """Read ``manifest.version`` from the pipeline's ``nextflow.config``.

    Pipeline mode gets the version from Nextflow itself (``--version``); this is
    the standalone-mode fallback, and keeps the two in step without a second
    place to bump.  ``nextflow.config`` sits one level above ``bin/``, where this
    script lives.  Returns ``None`` when it cannot be read.
    """
    config_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "nextflow.config"
    )
    try:
        with open(config_path, "r") as fh:
            content = fh.read()
    except Exception:
        return None
    match = re.search(
        r"manifest\s*\{.*?\bversion\s*=\s*['\"]([^'\"]+)['\"]",
        content, re.DOTALL,
    )
    return match.group(1) if match else None


def _meta_field(raw: Dict[str, Any], *keys: str) -> str:
    """First usable string among *keys* of a ``params_*.json`` dump.

    ``nextflow.config`` declares ``params.version = false`` to back Nextflow's
    own ``--version`` flag, and ``dumpParametersToJSON`` writes that boolean out
    as-is.  Passing it through ``str()`` is what made the report header read
    "vFalse", so booleans -- and their stringified spellings -- are skipped here
    and the caller falls back to ``manifest.version``.
    """
    for key in keys:
        value = raw.get(key)
        if value is None or isinstance(value, bool):
            continue
        text = str(value).strip()
        if text and text.lower() not in {"true", "false", "null", "none"}:
            return text
    return ""


# ---------------------------------------------------------------------------
# Curve thinning
# ---------------------------------------------------------------------------

def _log_spaced_indices(n: int, max_points: int, head: int = 0) -> List[int]:
    """Indices into a length-*n* sequence, thinned uniformly in ``log10(rank)``.

    The first *head* indices and the final index are always kept; the remaining
    budget is spread evenly in log10 of the rank, which is the axis these curves
    are drawn against.  The result is sorted and free of duplicates, so it may be
    shorter than *max_points*.

    Returns every index when the sequence already fits within the budget.
    """
    if n <= 0:
        return []
    if n <= max_points:
        return list(range(n))

    head = max(0, min(head, n))
    keep = set(range(head))
    keep.add(n - 1)

    budget = max(max_points - len(keep), 2)
    lo = math.log10(max(head, 1))
    hi = math.log10(n)
    for i in range(budget):
        rank = 10.0 ** (lo + (hi - lo) * i / (budget - 1))
        keep.add(min(max(int(round(rank)) - 1, 0), n - 1))
    return sorted(keep)


def build_knee_payload(path: Optional[str]) -> Optional[Dict[str, Any]]:
    """Read a ``UMIperCellSorted.txt`` into a thinned, plot-ready barcode-rank curve.

    The file holds one UMI count per barcode, which for a deeply sequenced sample
    is millions of lines.  Only ``_KNEE_MAX_POINTS`` of them are embedded, chosen
    by :func:`_log_spaced_indices`, and each retained point carries its true rank
    plus the cumulative UMI total up to that rank.  Carrying the exact rank and
    cumulative sum is what makes the thinning safe: the report reconstructs both
    the barcode-rank plot and the client-side second-derivative curve from these
    points without ever needing the ones in between.

    Returns ``None`` when the file is missing or holds no counts.  The payload:

    ``ranks`` / ``umis`` / ``cum``
        Parallel arrays over the retained barcodes (1-based rank, UMI count,
        cumulative UMIs up to and including that rank).
    ``n_barcodes``
        Total barcodes in the file, before thinning.
    ``n_above_min``
        Barcodes with at least ``min_umis`` UMIs -- the length of the curve the
        second-derivative cell calling operates on.
    ``min_umis``
        The floor that count was taken at.
    """
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path, "r") as fh:
            umis = [int(s) for s in (line.strip() for line in fh) if s.isdigit()]
    except Exception:
        return None
    if not umis:
        return None

    # The file is emitted sorted, but the rank axis only means anything if it is
    umis.sort(reverse=True)
    n = len(umis)
    n_above_min = sum(1 for v in umis if v >= _SD_MIN_UMIS)

    keep = set(_log_spaced_indices(n, _KNEE_MAX_POINTS, _KNEE_HEAD_RANKS))
    # The last barcode above the floor bounds the second-derivative curve, so the
    # client reaches the same end point as if it had the untruncated data
    if n_above_min:
        keep.add(n_above_min - 1)

    ranks: List[int] = []
    values: List[int] = []
    cum: List[int] = []
    running = 0
    for i, val in enumerate(umis):
        running += val
        if i in keep:
            ranks.append(i + 1)
            values.append(val)
            cum.append(running)

    return {
        "ranks":       ranks,
        "umis":        values,
        "cum":         cum,
        "n_barcodes":  n,
        "n_above_min": n_above_min,
        "min_umis":    _SD_MIN_UMIS,
    }


def thin_secondderiv_curve(entry: Dict[str, Any]) -> Dict[str, Any]:
    """Thin the cumulative-UMI curve of a ``*_knee_data.json`` payload.

    ``logX`` / ``logY`` / ``customdata`` carry one point per barcode above the
    UMI floor; the derivative arrays are already resampled onto a fixed grid by
    ``secondderiv_cellcalling.py`` and are left alone.  Points are picked in
    log10(rank), the curve's own x axis, and the floats are rounded to 6 decimals
    -- far below one screen pixel on a plot spanning a handful of log units.
    """
    log_x = entry.get("logX")
    if not isinstance(log_x, list) or not log_x:
        return entry

    log_y      = entry.get("logY") or []
    customdata = entry.get("customdata") or []
    idx = _log_spaced_indices(len(log_x), _SD_MAX_POINTS, _KNEE_HEAD_RANKS)

    thinned = dict(entry)
    thinned["logX"] = [round(log_x[i], 6) for i in idx]
    if len(log_y) == len(log_x):
        thinned["logY"] = [round(log_y[i], 6) for i in idx]
    if len(customdata) == len(log_x):
        thinned["customdata"] = [customdata[i] for i in idx]
    if isinstance(entry.get("derivX"), list):
        thinned["derivX"] = [round(v, 6) for v in entry["derivX"]]
    if isinstance(entry.get("derivY"), list):
        thinned["derivY"] = [round(v, 6) for v in entry["derivY"]]
    return thinned


def build_per_cell_payload(
    path: Optional[str],
    umi_threshold: Any,
) -> Tuple[Optional[Dict[str, List]], Optional[Dict[str, int]]]:
    """Load a per-cell metrics JSON and cut it down to what the scatters plot.

    The file is column-oriented with one entry per barcode carrying at least one
    read -- millions of them on a deeply sequenced sample, and by far the largest
    contributor to the report's size.  Two reductions are applied:

    * the ``Cell`` barcode strings are dropped, since no plot reads them;
    * every called cell is kept, while the non-cell background cloud is thinned
      by a fixed stride down to ``_PERCELL_MAX_NONCELLS`` points.

    The six scatters draw the non-cells as 4px translucent markers, so past a few
    tens of thousands they overplot into a solid shape whose outline and density
    gradient a uniform stride preserves exactly.  Barcodes appear in whitelist
    order, which is unrelated to any metric, so the stride is an unbiased sample.

    Cells are identified by ``IsCell`` -- membership of the filtered matrix --
    falling back to ``TotalReads > umi_threshold``, matching what the report does.

    Returns ``(columns, sampling)``, where *sampling* records the non-cell counts
    before and after thinning so the plot legends can say what is shown, or
    ``(None, None)`` when nothing usable was read.
    """
    if not path or not os.path.exists(path):
        return None, None
    try:
        with open(path, "r") as fh:
            raw = json.load(fh)
    except Exception:
        return None, None
    if not isinstance(raw, dict) or not raw:
        return None, None

    columns = {k: v for k, v in raw.items()
               if k in _PERCELL_COLUMNS and isinstance(v, list)}
    if not columns:
        return None, None
    n_rows = min(len(v) for v in columns.values())
    if n_rows == 0:
        return None, None

    is_cell_col = columns.get("IsCell")
    total_col   = columns.get("TotalReads")
    threshold   = safe_float(umi_threshold)
    if is_cell_col is not None:
        cell_flags = [bool(v) for v in is_cell_col[:n_rows]]
    elif total_col is not None and threshold is not None:
        cell_flags = [(safe_float(v) or 0.0) > threshold for v in total_col[:n_rows]]
    else:
        cell_flags = [True] * n_rows

    n_noncells = sum(1 for flag in cell_flags if not flag)
    stride = 1
    if n_noncells > _PERCELL_MAX_NONCELLS:
        stride = math.ceil(n_noncells / _PERCELL_MAX_NONCELLS)

    keep: List[int] = []
    seen_noncells = 0
    for i, flag in enumerate(cell_flags):
        if flag:
            keep.append(i)
        else:
            if seen_noncells % stride == 0:
                keep.append(i)
            seen_noncells += 1

    out: Dict[str, List] = {}
    for name, values in columns.items():
        if name in _PERCELL_PCT_COLUMNS:
            out[name] = [round(v, 3) if isinstance(v, float) else v
                         for v in (values[i] for i in keep)]
        else:
            out[name] = [values[i] for i in keep]

    sampling = {
        "noncells_total":  n_noncells,
        "noncells_shown":  len(keep) - (n_rows - n_noncells),
        "cells_shown":     n_rows - n_noncells,
    }
    return out, sampling


def _safe_read_json(path: Optional[str]) -> Dict[str, Any]:
    """Read a JSON object from *path*; return ``{}`` on any failure."""
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path, "r") as fh:
            obj = json.load(fh)
        return obj if isinstance(obj, dict) else {}
    except Exception:
        return {}


def parse_secondderiv_knee(path: Optional[str], sample_id: str) -> Optional[Dict[str, Any]]:
    """Read a ``*_knee_data.json`` produced by ``secondderiv_cellcalling.py``.

    The file is keyed by sample ID at the top level; the sole payload is
    returned when the key does not match the analytical ID exactly.
    """
    payload = _safe_read_json(path)
    if not payload:
        return None
    if sample_id in payload:
        entry = payload[sample_id]
    elif len(payload) == 1:
        entry = next(iter(payload.values()))
    else:
        return None
    return entry if isinstance(entry, dict) else None


def parse_cell_meta_tsv(path: Optional[str]) -> Dict[str, Any]:
    """Parse a per-cell TSV/CSV and compute per-column mean/median statistics."""
    if not path or not os.path.exists(path):
        return {}
    try:
        # Support both tab and comma-delimited files
        with open(path, "r", newline="") as fh:
            header_line = fh.readline()
            delimiter = "\t" if "\t" in header_line else ","
            fh.seek(0)
            rows = [r for r in csv.DictReader(fh, delimiter=delimiter) if r]

        if not rows:
            return {}

        numeric: Dict[str, List[float]] = {}
        for row in rows:
            for key, raw in row.items():
                if raw is None or str(raw).strip() == "":
                    continue
                try:
                    numeric.setdefault(key, []).append(float(raw))
                except ValueError:
                    pass

        def _median(xs: List[float]) -> float:
            ys = sorted(xs)
            mid = len(ys) // 2
            return ys[mid] if len(ys) % 2 == 1 else (ys[mid - 1] + ys[mid]) / 2.0

        result: Dict[str, Any] = {"n_rows": len(rows)}
        for key, values in numeric.items():
            if len(values) >= 2:
                result[f"{key}__mean"]   = sum(values) / len(values)
                result[f"{key}__median"] = _median(values)
        return result
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# Standalone-mode discovery helpers
# ---------------------------------------------------------------------------

def _extract_base_id(analytical_id: str) -> str:
    """Strip a known tool suffix from *analytical_id* to recover the base sample ID.

    Examples::

        "Nvec-1C_starsolo"         → "Nvec-1C"
        "Nvec-1C_geneext_starsolo" → "Nvec-1C"
        "BCA015_sample_alevinfry"  → "BCA015_sample"
    """
    for suffix in _ANALYTICAL_SUFFIXES:
        if analytical_id.endswith(suffix):
            return analytical_id[: -len(suffix)]
    return analytical_id


def _latest_glob(pattern: str) -> Optional[str]:
    """Return the most-recently modified file matching *pattern*, or ``None``."""
    matches = glob.glob(pattern)
    return max(matches, key=os.path.getmtime) if matches else None


def _probe(path: str) -> Optional[str]:
    """Return *path* if it exists on disk, otherwise ``None``."""
    return path if os.path.exists(path) else None


def _find_sankey(directory: str, analytical_id: str) -> Optional[str]:
    """Return this sample's Kraken sankey HTML from *directory*, or ``None``.

    PAVIAN names its output after the Kraken report it is handed -- here
    ``<id>.k2report`` -- and the exact spelling varies with the Pavian version
    (``<id>.sankey.html``, ``<id>.k2report.sankey.html``, ...), so the file is
    matched by sample ID rather than by one hard-coded name.  The remainder of
    the name has to start with a separator, otherwise sample ``S1`` would claim
    ``S10``'s report.
    """
    for path in sorted(glob.glob(os.path.join(directory, "*.sankey.html"))):
        name = os.path.basename(path)
        if not name.startswith(analytical_id):
            continue
        if name[len(analytical_id):][:1] in ("_", ".", "-"):
            return path
    return None


def _discover_cell_filtering(
    result_dir: str,
    analytical_id: str,
    files: Dict[str, str],
) -> None:
    """Populate cellsweep entries in *files* when outputs are present.

    Searches :data:`_CELLSWEEP_ROOTS` under *result_dir* for per-sample
    subdirectories, stopping at the first hit.
    """
    for cs_root in _CELLSWEEP_ROOTS:
        base = os.path.join(result_dir, cs_root, analytical_id)
        table = _probe(os.path.join(base, f"{analytical_id}_top_ambient_genes.csv"))
        if table:
            files["cs_table"] = table
            for key, suffix in [
                ("cs_plot_contrib", "_ambient_hat_histogram.png"),
                ("cs_plot_umap",    "_umap_comparison.png"),
            ]:
                candidate = _probe(os.path.join(base, f"{analytical_id}{suffix}"))
                if candidate:
                    files[key] = candidate
            break


def _discover_geneext(result_dir: str) -> Tuple[Optional[str], Optional[str]]:
    """Locate the GeneExt report and log under ``gene_ext/`` (standalone mode).

    Both are named after the annotation GeneExt wrote, which the pipeline calls
    ``geneext.gtf``; they are matched by suffix so a differently named output is
    still found.  Returns ``(report, log)``, either of which may be ``None``.
    """
    gene_ext = os.path.join(result_dir, "gene_ext")
    if not os.path.isdir(gene_ext):
        return None, None
    report = _latest_glob(os.path.join(gene_ext, "*.Report.html"))
    log    = _latest_glob(os.path.join(gene_ext, "*.GeneExt.log"))
    return report, log


def _discover_starsolo_samples(
    result_dir: str,
    active_samples: List[Dict],
    file_map: Dict[str, Dict],
) -> None:
    """Populate *active_samples* and *file_map* from ``mapping_STARsolo/``.

    Each immediate subdirectory of ``mapping_STARsolo/`` is treated as one
    analytical sample.  File paths are resolved deterministically from the
    known STARsolo output layout; no globbing is required.
    """
    starsolo_root = os.path.join(result_dir, "mapping_STARsolo")
    if not os.path.isdir(starsolo_root):
        return

    for sample_dir in sorted(os.listdir(starsolo_root)):
        sample_path = os.path.join(starsolo_root, sample_dir)
        if not os.path.isdir(sample_path):
            continue

        analytical_id = sample_dir
        base_id       = _extract_base_id(analytical_id)
        files: Dict[str, str] = {}

        # STAR alignment log files
        for key, suffix in [
            ("star_log",      "_Log.final.out"),
            ("star_full_log", "_Log.out"),
        ]:
            candidate = _probe(os.path.join(sample_path, f"{analytical_id}{suffix}"))
            if candidate:
                files[key] = candidate

        # Solo.out: prefer GeneFull_Ex50pAS over Gene
        solo_out = os.path.join(sample_path, f"{analytical_id}_Solo.out")
        for gene_dir in _SOLO_GENE_DIRS:
            gene_path = os.path.join(solo_out, gene_dir)
            if not os.path.isdir(gene_path):
                continue
            for key, fname in [
                ("star_summary", "Summary.csv"),
                ("cell_stats",   "CellReads.stats"),
                ("knee",         "UMIperCellSorted.txt"),
            ]:
                # setdefault: first (preferred) hit wins
                candidate = _probe(os.path.join(gene_path, fname))
                if candidate:
                    files.setdefault(key, candidate)

            # Second-derivative cell calling (only present for that cellfilter_method)
            sd_path = os.path.join(gene_path, _SECONDDERIV_DIR)
            for key, fname in [
                ("sd_knee",  f"{analytical_id}_knee_data.json"),
                ("sd_stats", f"{analytical_id}_secondderiv_statistics.json"),
            ]:
                candidate = _probe(os.path.join(sd_path, fname))
                if candidate:
                    files.setdefault(key, candidate)

        # rRNA / mtDNA metrics
        mt = _probe(os.path.join(
            result_dir, "rRNA_mtDNA", f"{analytical_id}_mt_rrna_metrics.txt"
        ))
        if mt:
            files["mt_rrna"] = mt

        # Sequencing saturation outputs
        sat_dir = os.path.join(result_dir, "saturation", analytical_id)
        for key, fname in [
            ("sat_log", f"{analytical_id}_saturation.log"),
            ("sat_img", f"{analytical_id}_saturation.png"),
            ("res_img", f"{analytical_id}_saturation_residuals.png"),
        ]:
            candidate = _probe(os.path.join(sat_dir, fname))
            if candidate:
                files[key] = candidate

        # Per-cell metrics JSON
        pc = _probe(os.path.join(
            result_dir, "summary_results", "per-cell_metrics",
            f"{analytical_id}_metrics.json",
        ))
        if pc:
            files["per_cell"] = pc

        # Kraken sankey (optional; scanned across common directory names)
        for sankey_root in _SANKEY_ROOTS:
            candidate = _find_sankey(os.path.join(result_dir, sankey_root), analytical_id)
            if candidate:
                files["sankey"] = candidate
                break

        _discover_cell_filtering(result_dir, analytical_id, files)

        active_samples.append({"id": analytical_id, "base": base_id, "source": "starsolo"})
        file_map[analytical_id] = files

        print(f"  [STARsolo] discovered {analytical_id}", file=sys.stderr)


def _discover_alevinfry_samples(
    result_dir: str,
    active_samples: List[Dict],
    file_map: Dict[str, Dict],
) -> None:
    """Populate *active_samples* and *file_map* from ``mapping_alevin/``.

    Expected layout::

        mapping_alevin/
          {analytical_id}/
            {analytical_id}_counts/
              quant.json                  → af_quant
            {analytical_id}_run/
              aux_info/
                meta_info.json            → af_meta

    Each immediate subdirectory of ``mapping_alevin/`` is one analytical sample.
    """
    alevin_root = os.path.join(result_dir, "mapping_alevin")
    if not os.path.isdir(alevin_root):
        return

    for sample_dir in sorted(os.listdir(alevin_root)):
        sample_path = os.path.join(alevin_root, sample_dir)
        if not os.path.isdir(sample_path):
            continue

        analytical_id = sample_dir
        base_id       = _extract_base_id(analytical_id)
        files: Dict[str, str] = {}

        # quant.json  →  af_quant
        quant = _probe(os.path.join(
            sample_path, f"{analytical_id}_counts", "quant.json"
        ))
        if quant:
            files["af_quant"] = quant

        # meta_info.json  →  af_meta
        meta = _probe(os.path.join(
            sample_path, f"{analytical_id}_run", "aux_info", "meta_info.json"
        ))
        if meta:
            files["af_meta"] = meta

        # quants_mat_cols.txt  →  af_cols (authoritative gene count)
        cols = _probe(os.path.join(
            sample_path, f"{analytical_id}_counts", "alevin", "quants_mat_cols.txt"
        ))
        if cols:
            files["af_cols"] = cols

        # Second-derivative cell calling
        sd_path = os.path.join(
            sample_path, f"{analytical_id}_counts", "alevin", _SECONDDERIV_DIR
        )
        for key, fname in [
            ("knee",     f"{analytical_id}_alevin_UMIperCellSorted.txt"),
            ("sd_knee",  f"{analytical_id}_knee_data.json"),
            ("sd_stats", f"{analytical_id}_secondderiv_statistics.json"),
        ]:
            candidate = _probe(os.path.join(sd_path, fname))
            if candidate:
                files[key] = candidate

        # Per-cell metrics JSON (same location as STARsolo)
        pc = _probe(os.path.join(
            result_dir, "summary_results", "per-cell_metrics",
            f"{analytical_id}_metrics.json",
        ))
        if pc:
            files["per_cell"] = pc

        # rRNA / mtDNA metrics (when produced for alevin samples)
        mt = _probe(os.path.join(
            result_dir, "rRNA_mtDNA", f"{analytical_id}_mt_rrna_metrics.txt"
        ))
        if mt:
            files["mt_rrna"] = mt

        _discover_cell_filtering(result_dir, analytical_id, files)

        active_samples.append({"id": analytical_id, "base": base_id, "source": "alevin"})
        file_map[analytical_id] = files

        print(f"  [alevin-fry] discovered {analytical_id}", file=sys.stderr)


def discover_result_dir(
    result_dir: str,
) -> Tuple[
    List[Dict],          # active_samples
    Dict[str, Dict],     # file_map
    Optional[str],       # samplesheet path
    Optional[str],       # run_config path
    Dict[str, str],      # pipeline_meta (version, commit)
]:
    """Auto-discover all analytical samples and their result files.

    Scans ``mapping_STARsolo/`` and ``mapping_alevin/`` under *result_dir*.
    Pipeline metadata (samplesheet, run config, version, commit) is resolved
    from the most-recently modified files inside ``pipeline_info/``.

    Parameters
    ----------
    result_dir:
        Absolute or relative path to the pipeline result directory.

    Returns
    -------
    active_samples:
        List of ``{"id": ..., "base": ..., "source": ...}`` dicts.
    file_map:
        ``{analytical_id: {file_key: path}}`` mapping.
    samplesheet_path:
        Path to the most-recent samplesheet CSV, or ``None``.
    run_config_path:
        Path to the most-recent run-config text file, or ``None``.
    pipeline_meta:
        Dict with ``"version"`` and ``"commit"`` extracted from the most-recent
        params JSON, or empty strings when unavailable.
    """
    active_samples: List[Dict]      = []
    file_map:       Dict[str, Dict] = {}

    # ── Pipeline info files ──────────────────────────────────────────────────
    pipeline_info    = os.path.join(result_dir, "pipeline_info")
    samplesheet_path = _latest_glob(os.path.join(pipeline_info, "samplesheet_*.csv"))
    run_config_path  = _latest_glob(os.path.join(pipeline_info, "run_config_*.txt"))
    latest_params    = _latest_glob(os.path.join(pipeline_info, "params_*.json"))

    pipeline_meta: Dict[str, str] = {}
    if latest_params:
        raw = _safe_read_json(latest_params)
        pipeline_meta["version"] = _meta_field(raw, "pipeline_version", "version")
        pipeline_meta["commit"]  = _meta_field(raw, "git_commit",      "commit")

    # ── Sample discovery ─────────────────────────────────────────────────────
    _discover_starsolo_samples(result_dir, active_samples, file_map)
    _discover_alevinfry_samples(result_dir, active_samples, file_map)

    return active_samples, file_map, samplesheet_path, run_config_path, pipeline_meta


# ---------------------------------------------------------------------------
# Pipeline-mode file mapping
# ---------------------------------------------------------------------------

def _build_file_map_from_cli(
    args: argparse.Namespace,
    active_samples: List[Dict],
) -> Dict[str, Dict]:
    """Map flat CLI file lists back to per-sample dicts (pipeline mode).

    Files are matched to sample IDs by checking whether the sample ID appears
    in the file's basename.  IDs are tested longest-first to avoid partial
    matches (e.g. ``sample1_subsampled`` before ``sample1``).
    """
    file_map: Dict[str, Dict]  = {s["id"]: {} for s in active_samples}
    source_by_id               = {s["id"]: s["source"] for s in active_samples}
    sorted_ids                 = sorted(file_map, key=len, reverse=True)

    def _map(cli_list: Optional[List[str]], key: str, allowed: Optional[set] = None) -> None:
        if not cli_list:
            return
        for path in cli_list:
            if not path:
                continue
            basename = os.path.basename(path)
            for sid in sorted_ids:
                if sid not in basename:
                    continue
                if allowed and source_by_id[sid] not in allowed:
                    continue
                file_map[sid][key] = path
                break

    # STARsolo-specific outputs
    _map(args.star_logs,               "star_log",      {"starsolo"})
    _map(args.star_summaries,          "star_summary",  {"starsolo"})
    _map(args.star_full_logs,          "star_full_log", {"starsolo"})
    _map(args.mt_rrna_metrics,         "mt_rrna",       {"starsolo"})
    _map(args.saturation_logs,         "sat_log",       {"starsolo"})
    _map(args.cell_stats,              "cell_stats",    {"starsolo"})
    _map(args.saturation_imgs,         "sat_img",       {"starsolo"})
    _map(args.residuals_imgs,          "res_img",       {"starsolo"})
    _map(args.sankey_files,            "sankey",        {"starsolo"})
    _map(args.per_cell_files,          "per_cell",      {"starsolo"})

    # Shared outputs by both mappers
    _map(args.knee_files,              "knee")
    _map(args.secondderiv_knee,        "sd_knee")
    _map(args.secondderiv_stats,       "sd_stats")

    # Alevin-specific output
    _map(args.af_meta_info,            "af_meta",       {"alevin"})
    _map(args.af_quant_json,           "af_quant",      {"alevin"})
    _map(args.af_cell_meta,            "af_cell",       {"alevin"})
    _map(args.af_mat_cols,             "af_cols",       {"alevin"})

    # Cell-filtering outputs (shared by both mappers)
    _map(args.cellsweep_tables,        "cs_table")
    _map(args.cellsweep_plots_contrib, "cs_plot_contrib")
    _map(args.cellsweep_plots_umap,    "cs_plot_umap")

    return file_map


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    # ── Mode selection & validation ──────────────────────────────────────────
    using_standalone = bool(args.result_dir)
    using_pipeline   = bool(args.analytical_samples)

    if not using_standalone and not using_pipeline:
        sys.exit(
            "Error: supply either --result-dir (standalone mode) "
            "or --analytical_samples + --samplesheet (pipeline mode)."
        )
    if using_standalone and using_pipeline:
        sys.exit(
            "Error: --result-dir and --analytical_samples are mutually exclusive."
        )

    # ── Initialise the three data sources that differ between modes ──────────
    if using_standalone:
        print(f"Standalone mode: scanning {args.result_dir}", file=sys.stderr)
        active_samples, file_map, samplesheet_path, run_config_path, pipeline_meta = \
            discover_result_dir(args.result_dir)

        if not active_samples:
            sys.exit(
                f"Error: no analytical samples discovered under {args.result_dir!r}. "
                "Check that mapping_STARsolo/ or mapping_alevin/ subdirectories exist."
            )

        raw_run_config = parse_run_config(run_config_path)
        version = args.version or pipeline_meta.get("version") or manifest_version() or "unknown"
        commit  = args.commit  or pipeline_meta.get("commit")  or "unknown"

        samplesheet_config: Dict[str, Any] = {}
        if samplesheet_path and os.path.exists(samplesheet_path):
            with open(samplesheet_path, "r") as fh:
                for row in csv.DictReader(fh):
                    if row.get("sample"):
                        samplesheet_config[row["sample"]] = row

        geneext_report, geneext_log = _discover_geneext(args.result_dir)

        print(
            f"Standalone mode: discovered {len(active_samples)} analytical samples.",
            file=sys.stderr,
        )

    else:  # pipeline mode
        raw_run_config = parse_run_config(args.run_config)
        version = args.version or manifest_version() or "unknown"
        commit  = args.commit  or "unknown"

        samplesheet_config = {}
        if args.samplesheet and os.path.exists(args.samplesheet):
            with open(args.samplesheet, "r") as fh:
                for row in csv.DictReader(fh):
                    if row.get("sample"):
                        samplesheet_config[row["sample"]] = row

        active_samples = []
        if os.path.exists(args.analytical_samples):
            with open(args.analytical_samples, "r") as fh:
                for row in csv.DictReader(fh):
                    active_samples.append({
                        "id":     row["analytical_id"],
                        "base":   row["base_id"],
                        "source": row.get("source", "auto").strip().lower(),
                    })
        else:
            sys.exit(f"Error: manifest not found: {args.analytical_samples}")

        # GeneExt runs once per pipeline run, on the merged alignments of every
        # sample, so its outputs are run-level and are not part of the file map
        geneext_report = (args.geneext_report or [None])[0]
        geneext_log    = (args.geneext_log    or [None])[0]

        file_map = _build_file_map_from_cli(args, active_samples)
        print(
            f"Pipeline mode: loaded {len(active_samples)} samples from manifest.",
            file=sys.stderr,
        )

    # ── Report metadata ──────────────────────────────────────────────────────
    report_metadata = {
        "project":  args.project,
        "pipeline": args.pipeline,
        "version":  version,
        "date":     datetime.date.today().isoformat(),
        "commit":   commit,
    }

    # ── Per-sample processing (identical for both modes) ─────────────────────
    # The overview mixes mappers in one table, so a single header cannot carry a
    # per-row definition: "% Mapped Reads" is uniquely-mapped for STARsolo and
    # all mapped fragments for alevin-fry. The Mapper column makes which one
    # applies explicit for every row. When every row turns out to be STARsolo the
    # column has only one meaning, and the header is sharpened to say so after
    # the loop below, once the mappers are known.
    _MAPPED_COL = 2
    global_cols = [
        "Sample", "Mapper", "% Mapped Reads", "N cells", "Saturation",
        "Reads Needed for Target Saturation", "Noise (% UMIs non-cell barcodes)",
        "Median Transcripts Per Cell", "% Intronic Reads", "% rRNA in Unique reads",
        "% rRNA in multimappers all pos", "% mtDNA in Unique reads",
        "% mtDNA in multimappers all pos",
    ]
    global_rows:        List[List]         = []
    samples_json_list:  List[Dict]         = []
    per_cell_data:      Dict[str, Any]     = {}
    saturation_images:  Dict[str, Any]     = {}
    knee_data:          Dict[str, Any]     = {}
    cell_filtering_data: Dict[str, Any]   = {}
    secondderiv_data:   Dict[str, Any]     = {}
    mappers_seen:       set                = set()

    def get_val(source: Dict, key: str, default: Any = "N/A") -> Any:
        return source.get(key, default) if source else default

    for sample_info in active_samples:
        s_id    = sample_info["id"]
        base_id = sample_info["base"]
        files   = file_map[s_id]

        # print(f"Processing {s_id}...", file=sys.stderr)

        source      = sample_info.get("source", "auto")
        using_star  = source == "starsolo" if source != "auto" \
                      else bool(files.get("star_log") and files.get("star_summary"))
        using_alevin = not using_star
        mappers_seen.add("starsolo" if using_star else "alevin")

        mt_stats          = parse_mt_rrna_metrics(files.get("mt_rrna"))
        reads_07_sat_val  = parse_saturation_log(files.get("sat_log")) if using_star else "N/A"
        intronic_pct      = to_pct(parse_starsolo_intronic(files.get("cell_stats"))) if using_star else "N/A"

        umi_threshold = 0
        if using_star:
            full_log = files.get("star_full_log")
            if not full_log and files.get("star_log"):
                guessed = files["star_log"].replace("Log.final.out", "Log.out")
                if os.path.exists(guessed):
                    full_log = guessed
            umi_threshold = extract_umi_cutoff(full_log)

        # STARsolo's own nUMImin, kept aside because umi_threshold is replaced by
        # the second-derivative cutoff below whenever the matrices were re-filtered
        # on one. The report shows both, so the two cell calls stay comparable.
        # 0 means the value was not in Log.out, which is nothing to report.
        starsolo_cutoff: Optional[int] = umi_threshold if using_star and umi_threshold else None

        if using_star:
            star_stats = parse_star_log(files.get("star_log"))
            solo_stats = parse_starsolo_summary(files.get("star_summary"))

            n_input_reads = get_val(star_stats, "Number of input reads")
            n_unique      = get_val(star_stats, "Uniquely mapped reads number")
            pct_unique    = get_val(star_stats, "Uniquely mapped reads %")
            pct_multi     = get_val(star_stats, "Number of reads mapped to multiple loci")
            if "%" not in str(pct_multi) and n_input_reads != "N/A" and pct_multi != "N/A":
                pct_multi = get_val(star_stats, "% of reads mapped to multiple loci", pct_multi)

            pct_multi_too_many = get_val(star_stats, "% of reads mapped to too many loci")
            pct_short          = get_val(star_stats, "% of reads unmapped: too short")
            pct_other          = get_val(star_stats, "% of reads unmapped: other")
            pct_r1_q30         = to_pct(get_val(solo_stats, "Q30 Bases in CB+UMI"))
            pct_r2_q30         = to_pct(get_val(solo_stats, "Q30 Bases in RNA read"))
            saturation         = to_pct(get_val(solo_stats, "Sequencing Saturation"))

            n_cells     = get_val(solo_stats, "Estimated Number of Cells",
                                  get_val(solo_stats, "Cells Detected"))
            mean_reads  = get_val(solo_stats, "Mean Reads per Cell")
            median_umis = get_val(solo_stats, "Median UMI per Cell")
            median_genes = get_val(solo_stats, "Median GeneFull_Ex50pAS per Cell")
            total_genes  = get_val(solo_stats, "Total GeneFull_Ex50pAS Detected")

            frac_in_cells = get_val(solo_stats, "Fraction of Unique Reads in Cells")
            noise_pct = "N/A"
            if frac_in_cells != "N/A":
                try:
                    noise_pct = f"{(1.0 - float(frac_in_cells)) * 100:.2f}%"
                except (ValueError, TypeError):
                    pass

        else:  # alevin-fry
            af_meta  = _safe_read_json(files.get("af_meta"))
            af_quant = _safe_read_json(files.get("af_quant"))
            af_cell  = parse_cell_meta_tsv(files.get("af_cell"))

            total_reads  = af_meta.get("total_reads", af_meta.get("num_processed", "N/A"))
            mapping_rate = af_meta.get("mapping_rate", "N/A")
            if mapping_rate == "N/A" and "percent_mapped" in af_meta:
                mapping_rate = float(af_meta["percent_mapped"]) / 100.0

            mapped_reads = af_meta.get("num_mapped", "N/A")
            if mapped_reads == "N/A" and total_reads != "N/A" and mapping_rate != "N/A":
                try:
                    mapped_reads = int(float(total_reads) * float(mapping_rate))
                except (ValueError, TypeError):
                    pass

            # salmon's num_mapped counts mapped *fragments*, including
            # multimappers (resolved downstream by cr-like-em) -- it is not
            # STARsolo's uniquely-mapped count. The dashboard labels these
            # fields per mapper so the two are not read as the same quantity.
            n_input_reads = total_reads
            n_unique      = mapped_reads
            pct_unique    = to_pct(mapping_rate)
            n_cells       = af_quant.get("num_quantified_cells",
                              af_quant.get("num_cells",
                              af_meta.get("final_num_cbs", "N/A")))
            total_genes   = count_alevin_genes(af_quant, files.get("af_cols"))

            mean_reads = "N/A"
            if total_reads != "N/A" and n_cells != "N/A":
                try:
                    if int(n_cells) > 0:
                        mean_reads = int(float(total_reads) / int(n_cells))
                except (ValueError, TypeError):
                    pass

            pct_multi = pct_multi_too_many = pct_short = pct_other = "N/A"
            pct_r1_q30 = pct_r2_q30 = saturation = noise_pct = "N/A"

            median_umis  = af_cell.get("nUMI__median",
                             af_cell.get("umis__median",
                             af_meta.get("mean_umis_per_cell", "N/A")))
            median_genes = af_cell.get("nGene__median",
                             af_cell.get("genes__median",
                             af_meta.get("mean_genes_per_cell", "N/A")))

        # ── Second-derivative cell calling ───────────────────────────────────
        # When cellfilter_method = "second_derivative", the cells are re-called on the
        # raw matrix: STARsolo's own filtered matrix is neither produced nor published,
        # and alevin-fry's quantified matrix is filtered again after its knee. Every
        # cell-level number below is recomputed on the matrix that was actually
        # filtered, so those values replace the ones the mapper reported for its own
        # cell set (STARsolo's Summary.csv, alevin-fry's quant.json / cell_meta.tsv).
        sd_cutoff: Optional[int] = None
        # alevin-fry's own gene count is the matrix width, i.e. the reference size;
        # re-calling cells yields a genes-detected count instead, which the dashboard
        # has to label differently.
        total_genes_scope = "reference" if using_alevin else "detected"
        sd_stats = _safe_read_json(files.get("sd_stats"))
        if sd_stats:
            sd_cutoff     = sd_stats.get("umi_threshold_applied")
            umi_threshold = sd_cutoff if sd_cutoff is not None else umi_threshold
            n_cells       = sd_stats.get("estimated_cells", n_cells)
            median_umis   = int(round(sd_stats["median_umis_per_cell"]))   if "median_umis_per_cell"  in sd_stats else median_umis
            median_genes  = int(round(sd_stats["median_genes_per_cell"]))  if "median_genes_per_cell" in sd_stats else median_genes
            if "total_genes_detected" in sd_stats:
                total_genes       = sd_stats["total_genes_detected"]
                total_genes_scope = "detected"

            # Reads per cell and noise both depend on the cell set, so recompute them
            n_reads = safe_float(n_input_reads)
            if n_reads and n_cells:
                try:
                    mean_reads = int(n_reads / int(n_cells))
                except (ValueError, TypeError, ZeroDivisionError):
                    pass

            frac_in_cells_sd = safe_float(sd_stats.get("fraction_unique_reads_in_cells"))
            if frac_in_cells_sd is not None:
                noise_pct = f"{(1.0 - frac_in_cells_sd) * 100:.2f}%"

            # Saturation is the duplication rate of the reads counted into the cells,
            # so it too describes the cell set rather than the library
            saturation_sd = safe_float(sd_stats.get("sequencing_saturation"))
            if saturation_sd is not None:
                saturation = to_pct(saturation_sd)

        sd_knee = parse_secondderiv_knee(files.get("sd_knee"), s_id)
        if sd_knee:
            secondderiv_data[s_id] = thin_secondderiv_curve(sd_knee)
            if sd_cutoff is None:
                sd_cutoff = sd_knee.get("threshold_umi")

        rrna_pct           = to_pct(get_val(mt_stats, "Percentage of rRNA reads (of uniquely mapped reads)"))
        rrna_multi_all     = to_pct(get_val(mt_stats, "Percentage of rRNA in multimapped reads (all alignments)"))
        rrna_multi_primary = to_pct(get_val(mt_stats, "Percentage of rRNA in multimapped reads (primary alignment)"))
        mtdna_unique       = to_pct(get_val(mt_stats, "Percentage of mtDNA reads (of mapped reads)"))
        mtdna_multi_all    = to_pct(get_val(mt_stats, "Percentage of mtDNA in multimapped reads (all alignments)"))
        mtdna_multi_primary = to_pct(get_val(mt_stats, "Percentage of mtDNA in multimapped reads (primary alignment)"))

        global_rows.append([
            s_id, "STARsolo" if using_star else "alevin-fry",
            pct_unique, fmt(n_cells), saturation, reads_07_sat_val, noise_pct,
            fmt(median_umis), intronic_pct, rrna_pct, rrna_multi_all,
            mtdna_unique, mtdna_multi_all,
        ])

        # ── Cell filtering ───────────────────────────────────────────────────
        cf_dict = None
        if files.get("cs_table"):
            cf_dict = {
                "method":    "cellsweep",
                "top_genes": parse_cell_filtering_table(files["cs_table"]),
                "plot1":     encode_image(files.get("cs_plot_contrib")),
                "plot2":     encode_image(files.get("cs_plot_umap")),
            }
        if cf_dict:
            cell_filtering_data[s_id] = cf_dict

        # ── Combined config (run config + samplesheet row for this sample) ───
        combined_config = raw_run_config.copy()
        if base_id in samplesheet_config:
            combined_config.update(samplesheet_config[base_id])

        # ── Per-cell metrics ─────────────────────────────────────────────────
        per_cell_cols, per_cell_sampling = build_per_cell_payload(
            files.get("per_cell"), umi_threshold
        )
        if per_cell_cols:
            per_cell_data[s_id] = per_cell_cols

        samples_json_list.append({
            "sample_id": s_id,
            "config":    combined_config,
            "source":    "starsolo" if using_star else "alevin",
            "mapping": {
                "n_reads_sample":       fmt(n_input_reads),
                "pct_r1_q30":           pct_r1_q30,
                "pct_r2_q30":           pct_r2_q30,
                "n_uniquely_mapped":    fmt(n_unique),
                "pct_uniquely_mapped":  pct_unique,
                "pct_multi_mapped":     pct_multi,
                "pct_multi_too_many":   pct_multi_too_many,
                "pct_unmapped_short":   pct_short,
                "pct_unmapped_other":   pct_other,
                "noise_pct":            noise_pct,
                "intronic_pct":         intronic_pct,
                "rrna_pct":             rrna_pct,
                "rrna_multi_all_pct":     rrna_multi_all,
                "rrna_multi_primary_pct": rrna_multi_primary,
                "mtdna_unique_pct":     mtdna_unique,
                "mtdna_multi_all_pct":  mtdna_multi_all,
                "mtdna_multi_primary_pct": mtdna_multi_primary,
                "n_cells":              fmt(n_cells),
                "saturation":           saturation,
                "reads_07_saturation":  reads_07_sat_val,
                "mean_reads_per_cell":  fmt(mean_reads),
                "median_umis":          fmt(median_umis),
                "median_genes":         fmt(median_genes),
                "total_genes":          fmt(total_genes),
                "total_genes_scope":    total_genes_scope,
            },
            "cell_calling": {
                "expected_cells": fmt(samplesheet_config.get(base_id, {}).get("expected_cells", "N/A")),
                "num_cells":      fmt(n_cells),
                "umi_threshold":  umi_threshold,
                "secondderiv_cutoff": sd_cutoff,
                "starsolo_cutoff":    starsolo_cutoff,
            },
            "taxonomy_sankey": extract_sankey_data(files.get("sankey")),
            "per_cell_sampling": per_cell_sampling,
        })

        # ── Knee data ────────────────────────────────────────────────────────
        knee_payload = build_knee_payload(files.get("knee"))
        if knee_payload:
            knee_data[s_id] = knee_payload

        # ── Saturation images ────────────────────────────────────────────────
        saturation_images[s_id] = {}
        if files.get("sat_img"):
            saturation_images[s_id]["saturation"] = encode_image(files["sat_img"])
        if files.get("res_img"):
            saturation_images[s_id]["residuals"]  = encode_image(files["res_img"])

    # Every row came from STARsolo, so the mapped-reads column carries a single
    # definition and the header can name it. A mixed run keeps the generic wording.
    if global_rows and mappers_seen == {"starsolo"}:
        global_cols[_MAPPED_COL] = "% Uniquely Mapped Reads"

    # ── Gene extension (run-level) ───────────────────────────────────────────
    geneext_data = build_geneext_payload(geneext_report, geneext_log)
    if geneext_data:
        print(
            f"GeneExt statistics read from the {geneext_data['source']}.",
            file=sys.stderr,
        )

    # ── Inject into HTML template ────────────────────────────────────────────
    # The template carries non-ASCII text (superscripts, arrows, primes), so the
    # encoding is pinned rather than taken from the locale -- a C/POSIX locale
    # would otherwise fail to read the template at all.
    with open(args.template, "r", encoding="utf-8") as fh:
        html_content = fh.read()

    # The payloads are parsed by the report, never read by a human, so they are
    # serialised without indentation: pretty-printing puts every element of the
    # per-cell and knee arrays on its own padded line, which roughly doubles the
    # size of the largest blocks for no benefit.
    def _dump(obj: Any) -> str:
        return json.dumps(obj, separators=(",", ":"))

    replacements = {
        "__RUN_METADATA_PLACEHOLDER__":      _dump(report_metadata),
        "__GLOBAL_DATA_PLACEHOLDER__":       _dump(
            {"overview": {"columns": global_cols, "rows": global_rows}}
        ),
        "__SAMPLES_DATA_PLACEHOLDER__":      _dump({"samples": samples_json_list}),
        "__PER_CELL_DATA_PLACEHOLDER__":     _dump(per_cell_data),
        "__SATURATION_IMAGES_PLACEHOLDER__": _dump(saturation_images),
        "__KNEE_DATA_PLACEHOLDER__":         _dump(knee_data),
        "__SECONDDERIV_DATA_PLACEHOLDER__":  _dump(secondderiv_data),
        "__CELLFILTERING_DATA_PLACEHOLDER__": _dump(cell_filtering_data),
        "__GENEEXT_DATA_PLACEHOLDER__":       _dump(geneext_data),
    }

    for placeholder, json_str in replacements.items():
        if placeholder in html_content:
            html_content = html_content.replace(placeholder, json_str)
        else:
            safe_ph = re.escape(placeholder)
            html_content = re.sub(r"\s*" + safe_ph + r"\s*", "\n" + json_str + "\n", html_content)

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(html_content)

    print(
        f"Successfully generated {args.output} "
        f"for {len(global_rows)} analytical samples.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
