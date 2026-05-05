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

# Candidate locations for optional cell-filtering outputs (checked in order).
_CELLSWEEP_ROOTS:  List[str] = ["cell_filtering/cellsweep", "cellsweep"]
_CELLBENDER_ROOTS: List[str] = ["cell_filtering/cellbender", "cellbender"]

# Candidate locations for Kraken sankey HTML files.
_SANKEY_ROOTS: List[str] = ["kraken", "taxonomy"]


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

    # ── Pipeline mode: visualisation file lists ──────────────────────────────
    parser.add_argument("--sankey_files",    nargs="*", help="*_kraken.sankey.html files")
    parser.add_argument("--per_cell_files",  nargs="*", help="Per-cell JSON files")
    parser.add_argument("--saturation_imgs", nargs="*", help="Saturation PNG images")
    parser.add_argument("--residuals_imgs",  nargs="*", help="Residuals PNG images")
    parser.add_argument("--saturation_logs",nargs="*", help="saturation.log files")
    parser.add_argument("--mt_rrna_metrics",nargs="*", help="*_mt_rrna_metrics.txt files")
    parser.add_argument("--knee_files",      nargs="*", help="UMIperCellSorted.txt files")

    parser.add_argument("--cellsweep_tables",        nargs="*")
    parser.add_argument("--cellsweep_plots_contrib", nargs="*")
    parser.add_argument("--cellsweep_plots_umap",    nargs="*")
    parser.add_argument("--cellbender_tables",        nargs="*")
    parser.add_argument("--cellbender_plots1",        nargs="*")
    parser.add_argument("--cellbender_plots2",        nargs="*")

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
    """Extract the htmlwidgets JSON payload from a sankey HTML file."""
    if not html_path or not os.path.exists(html_path):
        return None
    try:
        with open(html_path, "r") as fh:
            content = fh.read()
        match = re.search(
            r'<script type="application/json" data-for="htmlwidget-[a-f0-9]+">(.*?)</script>',
            content, re.DOTALL,
        )
        if match:
            return json.loads(match.group(1)).get("x")
    except Exception:
        pass
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


def parse_cell_filtering_table(path: Optional[str], method: str) -> List[Dict[str, Any]]:
    """Parse a cellsweep or cellbender output table into a list of row dicts."""
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
                if method == "cellsweep":
                    data.append({
                        "gene":        cleaned.get("gene", "N/A"),
                        "ambient_hat": safe_float(cleaned.get("ambient_hat")),
                    })
                elif method == "cellbender":
                    data.append({
                        "gene_name":        cleaned.get("gene_name", "N/A"),
                        "n_removed":        safe_float(cleaned.get("n_removed")),
                        "fraction_removed": safe_float(cleaned.get("fraction_removed")),
                    })
        sort_key = "ambient_hat" if method == "cellsweep" else "n_removed"
        data.sort(key=lambda x: x[sort_key] if x[sort_key] is not None else -1, reverse=True)
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

    Returns the fraction as a ``float``, or ``"N/A"`` when unavailable.
    """
    if not path or not os.path.exists(path):
        return "N/A"
    try:
        mapped_reads = 0
        intronic_sum = 0
        with open(path, "r") as fh:
            next(fh)  # skip header
            for line in fh:
                parts = line.strip().split()
                if len(parts) < 6:
                    continue
                intronic_sum += int(parts[4]) + int(parts[5])
                mapped_reads += int(parts[2]) + int(parts[3]) + int(parts[4]) + int(parts[5])
        if mapped_reads > 0:
            return intronic_sum / mapped_reads
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


def _discover_cell_filtering(
    result_dir: str,
    analytical_id: str,
    files: Dict[str, str],
) -> None:
    """Populate cellsweep / cellbender entries in *files* when outputs are present.

    Searches :data:`_CELLSWEEP_ROOTS` and :data:`_CELLBENDER_ROOTS` under
    *result_dir* for per-sample subdirectories, stopping at the first hit for
    each method.
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

    for cb_root in _CELLBENDER_ROOTS:
        base = os.path.join(result_dir, cb_root, analytical_id)
        table = _probe(os.path.join(base, f"{analytical_id}_cellbender.tsv"))
        if table:
            files["cb_table"] = table
            for key, suffix in [
                ("cb_plot1", "_plot1.png"),
                ("cb_plot2", "_plot2.png"),
            ]:
                candidate = _probe(os.path.join(base, f"{analytical_id}{suffix}"))
                if candidate:
                    files[key] = candidate
            break


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
            candidate = _probe(os.path.join(
                result_dir, sankey_root, f"{analytical_id}_kraken.sankey.html"
            ))
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
        pipeline_meta["version"] = str(raw.get("pipeline_version", raw.get("version", "")))
        pipeline_meta["commit"]  = str(raw.get("git_commit",       raw.get("commit",  "")))

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
    _map(args.knee_files,              "knee",          {"starsolo"})

    _map(args.af_meta_info,            "af_meta",       {"alevin"})
    _map(args.af_quant_json,           "af_quant",      {"alevin"})
    _map(args.af_cell_meta,            "af_cell",       {"alevin"})

    # Cell-filtering outputs are shared (either tool can emit them)
    _map(args.cellsweep_tables,        "cs_table")
    _map(args.cellsweep_plots_contrib, "cs_plot_contrib")
    _map(args.cellsweep_plots_umap,    "cs_plot_umap")
    _map(args.cellbender_tables,       "cb_table")
    _map(args.cellbender_plots1,       "cb_plot1")
    _map(args.cellbender_plots2,       "cb_plot2")

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
        version = args.version or pipeline_meta.get("version") or "unknown"
        commit  = args.commit  or pipeline_meta.get("commit")  or "unknown"

        samplesheet_config: Dict[str, Any] = {}
        if samplesheet_path and os.path.exists(samplesheet_path):
            with open(samplesheet_path, "r") as fh:
                for row in csv.DictReader(fh):
                    if row.get("sample"):
                        samplesheet_config[row["sample"]] = row

        print(
            f"Standalone mode: discovered {len(active_samples)} analytical samples.",
            file=sys.stderr,
        )

    else:  # pipeline mode
        raw_run_config = parse_run_config(args.run_config)
        version = args.version or "unknown"
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
    global_cols = [
        "Sample", "% Uniquely Mapped Reads", "N cells", "Saturation",
        "Reads Needed for Target Saturation", "Noise (% UMIs non-cell barcodes)",
        "Median Transcripts Per Cell", "% Intronic Reads", "% rRNA in Unique reads",
        "% mtDNA in Unique reads", "% mtDNA in multimappers all pos",
    ]
    global_rows:        List[List]         = []
    samples_json_list:  List[Dict]         = []
    per_cell_data:      Dict[str, Any]     = {}
    saturation_images:  Dict[str, Any]     = {}
    knee_data:          Dict[str, Any]     = {}
    cell_filtering_data: Dict[str, Any]   = {}

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

            n_input_reads = total_reads
            n_unique      = mapped_reads
            pct_unique    = to_pct(mapping_rate)
            n_cells       = af_quant.get("num_quantified_cells",
                              af_quant.get("num_cells",
                              af_meta.get("final_num_cbs", "N/A")))
            total_genes   = af_quant.get("num_genes",
                              af_quant.get("num_targets", "N/A"))

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

        rrna_pct           = to_pct(get_val(mt_stats, "Percentage of rRNA reads (of uniquely mapped reads)"))
        mtdna_unique       = to_pct(get_val(mt_stats, "Percentage of mtDNA reads (of mapped reads)"))
        mtdna_multi_all    = to_pct(get_val(mt_stats, "Percentage of mtDNA in multimapped reads (all alignments)"))
        mtdna_multi_primary = to_pct(get_val(mt_stats, "Percentage of mtDNA in multimapped reads (primary alignment)"))

        global_rows.append([
            s_id, pct_unique, fmt(n_cells), saturation, reads_07_sat_val, noise_pct,
            fmt(median_umis), intronic_pct, rrna_pct, mtdna_unique, mtdna_multi_all,
        ])

        # ── Cell filtering ───────────────────────────────────────────────────
        cf_dict = None
        if files.get("cs_table"):
            cf_dict = {
                "method":    "cellsweep",
                "top_genes": parse_cell_filtering_table(files["cs_table"], "cellsweep"),
                "plot1":     encode_image(files.get("cs_plot_contrib")),
                "plot2":     encode_image(files.get("cs_plot_umap")),
            }
        elif files.get("cb_table"):
            cf_dict = {
                "method":    "cellbender",
                "top_genes": parse_cell_filtering_table(files["cb_table"], "cellbender"),
                "plot1":     encode_image(files.get("cb_plot1")),
                "plot2":     encode_image(files.get("cb_plot2")),
            }
        if cf_dict:
            cell_filtering_data[s_id] = cf_dict

        # ── Combined config (run config + samplesheet row for this sample) ───
        combined_config = raw_run_config.copy()
        if base_id in samplesheet_config:
            combined_config.update(samplesheet_config[base_id])

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
            },
            "cell_calling": {
                "expected_cells": fmt(samplesheet_config.get(base_id, {}).get("expected_cells", "N/A")),
                "num_cells":      fmt(n_cells),
                "umi_threshold":  umi_threshold,
            },
            "taxonomy_sankey": extract_sankey_data(files.get("sankey")),
        })

        # ── Per-cell JSON ────────────────────────────────────────────────────
        if files.get("per_cell"):
            try:
                with open(files["per_cell"], "r") as fh:
                    per_cell_data[s_id] = json.load(fh)
            except Exception:
                pass

        # ── Knee data ────────────────────────────────────────────────────────
        knee_data[s_id] = []
        if files.get("knee"):
            try:
                with open(files["knee"], "r") as fh:
                    knee_data[s_id] = [
                        int(line.strip()) for line in fh if line.strip().isdigit()
                    ]
            except Exception:
                pass

        # ── Saturation images ────────────────────────────────────────────────
        saturation_images[s_id] = {}
        if files.get("sat_img"):
            saturation_images[s_id]["saturation"] = encode_image(files["sat_img"])
        if files.get("res_img"):
            saturation_images[s_id]["residuals"]  = encode_image(files["res_img"])

    # ── Inject into HTML template ────────────────────────────────────────────
    with open(args.template, "r") as fh:
        html_content = fh.read()

    replacements = {
        "__RUN_METADATA_PLACEHOLDER__":      json.dumps(report_metadata,      indent=2),
        "__GLOBAL_DATA_PLACEHOLDER__":       json.dumps(
            {"overview": {"columns": global_cols, "rows": global_rows}}, indent=2
        ),
        "__SAMPLES_DATA_PLACEHOLDER__":      json.dumps({"samples": samples_json_list}, indent=2),
        "__PER_CELL_DATA_PLACEHOLDER__":     json.dumps(per_cell_data,         indent=2),
        "__SATURATION_IMAGES_PLACEHOLDER__": json.dumps(saturation_images,     indent=2),
        "__KNEE_DATA_PLACEHOLDER__":         json.dumps(knee_data,             indent=2),
        "__CELLFILTERING_DATA_PLACEHOLDER__": json.dumps(cell_filtering_data,  indent=2),
    }

    for placeholder, json_str in replacements.items():
        if placeholder in html_content:
            html_content = html_content.replace(placeholder, json_str)
        else:
            safe_ph = re.escape(placeholder)
            html_content = re.sub(r"\s*" + safe_ph + r"\s*", "\n" + json_str + "\n", html_content)

    with open(args.output, "w") as fh:
        fh.write(html_content)

    print(
        f"Successfully generated {args.output} "
        f"for {len(global_rows)} analytical samples.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
