#!/usr/bin/env python
"""
Aggregate QC log metrics from STARsolo, ParseBio, CellRanger, sci-rocket and
alevin-fry pipelines into a single TSV with universal column names.

Optionally writes a JSON file with CellReads-derived percentages per sample (STARsolo-only).

Example
-------
Run on a results directory that contains subfolders like:

    mapping_STARsolo/
    ParseBio_pipeline/
    CellRanger_pipeline/
    sci-rocket/
    mapping_alevin/

    python dashboard_mappingstats.py /analysis/ \
        --outfile mapping_stats.tsv \
        --json cellreads_summary.json
"""


from __future__ import annotations

from pathlib import Path
import argparse
import json
import re
import subprocess
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_command_line_arguments(
    argv: Optional[List[str]] = None,
) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Aggregate mapping QC log_metrics into one TSV file."
    )

    parser.add_argument(
        "outdir",
        type=Path,
        help=(
            "Results directory containing mapping_STARsolo/, ParseBio_pipeline/, "
            "CellRanger_pipeline/, sci-rocket/, mapping_alevin/, …"
        ),
    )
    parser.add_argument(
        "-o",
        "--outfile",
        default=Path("mapping_stats.tsv"),
        type=Path,
        help="Destination TSV (default: mapping_stats.tsv)",
    )
    parser.add_argument(
        "--json",
        type=Path,
        help="Optional JSON with CellReads-derived percentages (STARsolo only).",
    )

    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def run_command(cmd: List[str]) -> str:
    """Run a shell command and return stdout, or '' on failure."""
    try:
        return subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"[WARNING] Command failed: {' '.join(cmd)} -> {e}")
        return ""


def clean_int(val: object) -> Optional[int]:
    """Parse an integer from a string that may contain commas or '%'."""
    try:
        return int(str(val).replace(",", "").replace("%", ""))
    except (TypeError, ValueError):
        return None


def clean_float(val: object) -> Optional[float]:
    """Parse a float from a string that may contain '%'."""
    try:
        return float(str(val).replace("%", ""))
    except (TypeError, ValueError):
        return None


def convert_to_pct(value: Optional[float]) -> str:
    """Format a fraction (0..1) as a percentage string with two decimals."""
    if value is None:
        return ""
    return f"{value * 100:.2f}%"


def safe_fraction(numerator: Optional[float], denominator: Optional[float]) -> Optional[float]:
    """Return numerator/denominator, safely."""
    if numerator is None or denominator in (None, 0):
        return None
    return float(numerator) / float(denominator)


def first_existing(paths: Iterable[Path]) -> Optional[Path]:
    """Return the first existing path in an iterable, or None."""
    for p in paths:
        if p.exists():
            return p
    return None


def count_file_lines(path: Path) -> Optional[int]:
    """Count number of lines in a text file, or None on failure."""
    if not path.exists():
        print(f"[WARNING] File not found (cannot count lines): {path}")
        return None
    try:
        with path.open("rb") as fh:
            return sum(1 for _ in fh)
    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to count lines in {path}: {e}")
        return None


# ---------------------------------------------------------------------------
# STARsolo helpers
# ---------------------------------------------------------------------------

_STAR_PATTERNS: Dict[str, Tuple[str, Optional[Callable[[str], object]]]] = {
    "N reads/sample": ("Number of input reads", clean_int),
    "N uniquely mapped reads": ("Uniquely mapped reads number", clean_int),
    "% uniquely mapped reads": ("Uniquely mapped reads %", None),
    "% multi-mapped reads": ("% of reads mapped to multiple loci", None),
    "% multi-mapped reads: too many": ("% of reads mapped to too many loci", None),
    "% unmapped: too short": ("% of reads unmapped: too short", None),
    "% unmapped: other": ("% of reads unmapped: other", None),
}


def parse_star_logfinal(log_file: Path) -> Dict[str, object]:
    """Extract selected metrics from STARsolo Log.final.out."""
    metrics: Dict[str, object] = {}

    if not log_file.exists():
        print(f"[WARNING] STARsolo Log.final.out not found: {log_file}")
        return metrics

    try:
        for line in log_file.read_text().splitlines():
            for col_name, (pattern, caster) in _STAR_PATTERNS.items():
                if pattern in line:
                    token = line.split()[-1]
                    metrics[col_name] = caster(token) if caster else token
    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to parse {log_file}: {e}")

    return metrics


def choose_gene_directory(solo_out: Path) -> Tuple[Path, str]:
    """
    Select STARsolo gene directory and its config label.
    Prefers GeneFull_Ex50pAS when present; otherwise returns Gene.
    """
    gene_full = solo_out / "GeneFull_Ex50pAS"
    if gene_full.exists():
        return gene_full, "GeneFull_Ex50pAS"
    return solo_out / "Gene", "Gene"


def load_summary_csv(csv_path: Path) -> pd.DataFrame:
    """Load STARsolo Summary.csv as a 2-column table with the first column as index."""
    if not csv_path.exists():
        print(f"[WARNING] Summary CSV file not found: {csv_path}")
        return pd.DataFrame()

    try:
        df = pd.read_csv(csv_path)
        df.set_index(df.columns[0], inplace=True)
        return df
    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Could not read STARsolo Summary.csv at {csv_path}: {e}")
        return pd.DataFrame()


def extract_star_summary_fields(df: pd.DataFrame, cfg_name: str) -> Dict[str, str]:
    """
    Pick and normalize selected fields from STARsolo Summary.csv.

    Metrics that represent proportions (Q30, Saturation) are formatted as
    percentage strings.
    """
    wanted = {
        "Q30 Bases in CB+UMI": "Q30 Bases in CB+UMI",
        "Q30 Bases in RNA read": "Q30 Bases in RNA read",
        "Saturation": "Sequencing Saturation",
        "Mean Reads per Cell": "Mean Reads per Cell",
        "Median UMI Counts per Cell": "Median UMI per Cell",
        "Median Genes per Cell": f"Median {cfg_name} per Cell",
        "Total Genes Detected": f"Total {cfg_name} Detected",
    }
    out: Dict[str, str] = {}
    if df.empty:
        return out

    for key, value in wanted.items():
        if value not in df.index:
            continue
        try:
            val = df.at[value, df.columns[-1]]
            if "Q30" in key or "Saturation" in key:
                out[key] = convert_to_pct(val)
            else:
                out[key] = str(int(val))
        except Exception:  # noqa: BLE001
            out[key] = ""

    return out


def extract_umi_cutoff(log_out: Path, cfg_name: str) -> Optional[int]:
    """Extract nUMImin (cell-calling cutoff) from STARsolo Log.out."""
    if not log_out.exists():
        print(f"[WARNING] STARsolo Log.out not found: {log_out}")
        return None

    try:
        in_block = False
        for line in log_out.read_text().splitlines():
            if f"Starting Solo post-map for {cfg_name}" in line:
                in_block = True
                continue
            if in_block and "Starting Solo post-map for" in line:
                break
            if in_block and "cellFiltering" in line:
                match = re.search(r"nUMImin=(\d+)", line)
                if match:
                    return int(match.group(1))
    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to parse {log_out}: {e}")
    return None


def find_star_file(sample_dir: Path, suffix: str) -> Optional[Path]:
    """
    Find a STARsolo file (Log.final.out, Log.out, Solo.out directory) that may
    have a sample prefix (e.g. r12345_Log.final.out).
    """
    direct = sample_dir / suffix
    if direct.exists():
        return direct

    matches = list(sample_dir.glob(f"*_{suffix}"))
    if matches:
        return matches[0]
    return None


# ---------------------------------------------------------------------------
# CellReads helpers (STARsolo)
# ---------------------------------------------------------------------------

def add_cellreads_metrics(read_stats: Path, barcodes: Path) -> Dict[str, str]:
    """
    Compute CellReads-derived percentages as strings with '%' suffix.

    Returns keys like:
        - pct_noise
        - pct_exonic_reads
        - pct_intronic_reads
        - pct_intergenic_reads
        - pct_mitochondrial_reads
        - pct_exonicAS_reads
        - pct_intronicAS_reads
    """
    if not read_stats.exists() or not barcodes.exists():
        print(f"[WARNING] CellReads stats or barcodes file not found: {read_stats}, {barcodes}")
        return {}

    try:
        stats_df = pd.read_csv(read_stats, sep="\t", index_col=0)
        filtered_cbs = set(pd.read_csv(barcodes, header=None)[0])

        total_cbmatch = stats_df.get("cbMatch", pd.Series(dtype=float)).sum()
        cbmatch_noise = (
            stats_df.loc["CBnotInPasslist", "cbMatch"]
            if "CBnotInPasslist" in stats_df.index
            else 0
        )

        filtered_df = stats_df.loc[stats_df.index.intersection(filtered_cbs)]
        summed = filtered_df.sum(numeric_only=True)
        genome_total = summed.get("genomeU", 0) + summed.get("genomeM", 0) or None

        out: Dict[str, str] = {}

        if total_cbmatch:
            out["pct_noise"] = convert_to_pct(cbmatch_noise / total_cbmatch)

        if genome_total:
            out["pct_exonic_reads"] = convert_to_pct(summed.get("exonic", 0) / genome_total)
            out["pct_intronic_reads"] = convert_to_pct(summed.get("intronic", 0) / genome_total)

            intergenic = (
                genome_total
                - summed.get("exonic", 0)
                - summed.get("intronic", 0)
                - summed.get("exonicAS", 0)
                - summed.get("intronicAS", 0)
            )
            out["pct_intergenic_reads"] = convert_to_pct(intergenic / genome_total)
            out["pct_mitochondrial_reads"] = convert_to_pct(summed.get("mito", 0) / genome_total)
            out["pct_exonicAS_reads"] = convert_to_pct(
                summed.get("exonicAS", 0) / genome_total
            )
            out["pct_intronicAS_reads"] = convert_to_pct(
                summed.get("intronicAS", 0) / genome_total
            )

        print("Extracted CellReads metrics:", out)
        return out

    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to compute CellReads metrics from {read_stats}: {e}")
        return {}


# ---------------------------------------------------------------------------
# sci-rocket helpers
# ---------------------------------------------------------------------------

def load_scirocket_qc(js_path: Path) -> Optional[Dict]:
    """
    Load qc_data.js from sci-rocket as a JSON-like dict.

    Expects:
        var data = { ... };
    """
    if not js_path.exists():
        print(f"[WARNING] sci-rocket qc_data.js not found: {js_path}")
        return None

    try:
        text = js_path.read_text()
        match = re.search(r"var\s+data\s*=\s*(\{.*\})\s*;?\s*$", text, re.DOTALL)
        if not match:
            print(f"[WARNING] Could not locate 'var data = {{...}}' block in {js_path}")
            return None
        json_str = match.group(1)
        return json.loads(json_str)
    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to parse sci-rocket qc_data.js at {js_path}: {e}")
        return None


def parse_scirocket_run(js_path: Path) -> List[Dict[str, object]]:
    """
    Parse one sci-rocket qc_data.js file and return a list of per-sample rows.
    """
    data = load_scirocket_qc(js_path)
    if not data:
        return []

    rows: List[Dict[str, object]] = []
    sample_succes = data.get("sample_succes", {})  # sic

    for sample_id, stats in sample_succes.items():
        row: Dict[str, object] = {
            "Sample": sample_id,
            "Software": "sci-rocket",
        }

        total_reads = stats.get("total_reads")
        if total_reads is not None:
            row["N reads/sample"] = int(total_reads)

        sat = stats.get("sequencing_saturation")
        if sat is not None:
            row["Saturation"] = convert_to_pct(float(sat))

        est_cells = stats.get("estimated_cells")
        if est_cells is not None:
            row["N cells"] = int(est_cells)

        mean_reads = stats.get("mean_reads_per_cell")
        if mean_reads is not None:
            row["Mean Reads per Cell"] = int(mean_reads)

        uniq_frac = stats.get("perc_unique_reads_genome_unique")
        if uniq_frac is not None:
            row["% uniquely mapped reads"] = convert_to_pct(float(uniq_frac))

        exonic = stats.get("total_exonic_reads", 0)
        intronic = stats.get("total_intronic_reads", 0)
        intergenic = stats.get("total_intergenic_reads", 0)
        mito = stats.get("total_mitochondrial_reads", 0)
        exonic_as = stats.get("total_exonicAS_reads", 0)
        intronic_as = stats.get("total_intronicAS_reads", 0)

        total_genome = exonic + intronic + intergenic + mito + exonic_as + intronic_as

        if total_genome > 0:
            row["% exonic reads"] = convert_to_pct(exonic / total_genome)
            row["% Intronic reads"] = convert_to_pct(intronic / total_genome)
            row["% intergenic reads"] = convert_to_pct(intergenic / total_genome)
            row["% mtDNA in Unique reads"] = convert_to_pct(mito / total_genome)
            row["% exonicAS reads"] = convert_to_pct(exonic_as / total_genome)
            row["% intronicAS reads"] = convert_to_pct(intronic_as / total_genome)

        rows.append(row)

    return rows


# ---------------------------------------------------------------------------
# alevin-fry helpers
# ---------------------------------------------------------------------------

def pick_column(df: pd.DataFrame, candidates: List[str]) -> Optional[pd.Series]:
    """Return the first existing column from candidates in df, or None."""
    for name in candidates:
        if name in df.columns:
            return df[name]
    return None


def parse_salmon_meta(meta_path: Path) -> Dict[str, Optional[int]]:
    """
    Parse Salmon/alevin-fry meta_info.json and return total/mapped read counts.
    """
    result: Dict[str, Optional[int]] = {
        "total_reads": None,
        "mapped_reads": None,
    }

    if not meta_path or not meta_path.exists():
        print(f"[WARNING] alevin-fry meta_info.json not found: {meta_path}")
        return result

    try:
        meta = json.loads(meta_path.read_text())
        total = (
            meta.get("num_processed")
            or meta.get("num_fragments")
            or meta.get("num_reads")
        )
        mapped = (
            meta.get("num_mapped")
            or meta.get("num_mapped_fragments")
            or meta.get("mapped")
        )
        result["total_reads"] = int(total) if total is not None else None
        result["mapped_reads"] = int(mapped) if mapped is not None else None
    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to parse alevin-fry meta_info.json at {meta_path}: {e}")

    return result


def parse_alevin_quant(quant_path: Path) -> Dict[str, Optional[int]]:
    """
    Parse alevin-fry quant.json and return:
        - n_cells: num_quantified_cells
        - num_genes: num_genes
    """
    summary: Dict[str, Optional[int]] = {
        "n_cells": None,
        "num_genes": None,
    }

    if not quant_path or not quant_path.exists():
        print(f"[WARNING] alevin-fry quant.json not found: {quant_path}")
        return summary

    try:
        data = json.loads(quant_path.read_text())
        if "num_quantified_cells" in data:
            summary["n_cells"] = int(data["num_quantified_cells"])
        if "num_genes" in data:
            summary["num_genes"] = int(data["num_genes"])
    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to parse alevin-fry quant.json at {quant_path}: {e}")

    return summary


def summarize_cell_meta(cell_meta_path: Path) -> Dict[str, Optional[int]]:
    """
    Summarize alevin-fry cell_meta.tsv-like file.

    Returns:
        - n_cells
        - mean_reads_per_cell
        - median_umi_per_cell
        - median_genes_per_cell
    """
    summary: Dict[str, Optional[int]] = {
        "n_cells": None,
        "mean_reads_per_cell": None,
        "median_umi_per_cell": None,
        "median_genes_per_cell": None,
    }

    if not cell_meta_path or not cell_meta_path.exists():
        print(f"[WARNING] alevin-fry cell_meta file not found: {cell_meta_path}")
        return summary

    try:
        df = pd.read_csv(cell_meta_path, sep="\t", comment="#")
    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to read alevin-fry cell_meta at {cell_meta_path}: {e}")
        return summary

    if df.empty:
        print(f"[WARNING] alevin-fry cell_meta file is empty: {cell_meta_path}")
        return summary

    summary["n_cells"] = int(df.shape[0])

    reads_series = pick_column(
        df,
        ["n_frags", "n_fragments", "n_frag", "nCount_RNA", "n_reads", "reads"],
    )
    if reads_series is not None:
        summary["mean_reads_per_cell"] = int(round(reads_series.mean()))

    umi_series = pick_column(df, ["nUMI", "n_umi", "umis", "UMI", "nCount_UMI"])
    if umi_series is not None:
        summary["median_umi_per_cell"] = int(round(umi_series.median()))

    gene_series = pick_column(df, ["nGene", "n_genes", "genes", "nGenes"])
    if gene_series is not None:
        summary["median_genes_per_cell"] = int(round(gene_series.median()))

    return summary


def parse_alevinfry_sample(sample_root: Path) -> Dict[str, object]:
    """
    Parse one alevin-fry sample directory under mapping_alevin/<sample_id>.

    Uses:
      - <sample_id>_run/aux_info/meta_info.json   → total & mapped reads
      - <sample_id>_counts/quant.json             → N cells, Total Genes Detected
      - cell_meta.tsv                             → per-cell summaries
    """
    sample_name = sample_root.name
    row: Dict[str, object] = {
        "Sample": sample_name,
        "Software": "alevin-fry",
        "N reads/sample": None,
        "N uniquely mapped reads": None,
        "% uniquely mapped reads": None,
        "N cells": None,
        "Mean Reads per Cell": None,
        "Median UMI Counts per Cell": None,
        "Median Genes per Cell": None,
        "Total Genes Detected": None,
    }

    # Identify *_run and *_counts directories
    run_dir = first_existing(d for d in sample_root.glob("*_run") if d.is_dir())
    counts_dir = first_existing(d for d in sample_root.glob("*_counts") if d.is_dir())

    # ---- meta_info.json → total & mapped reads ----
    meta_candidates: List[Path] = []
    if run_dir is not None:
        meta_candidates.extend(
            [
                run_dir / "aux_info" / "meta_info.json",  # your real path
                run_dir / "meta_info.json",
                run_dir / "quant" / "meta_info.json",
            ]
        )
    if counts_dir is not None:
        meta_candidates.extend(
            [
                counts_dir / "meta_info.json",
                counts_dir / "quant" / "meta_info.json",
            ]
        )
    meta_candidates.append(sample_root / "meta_info.json")
    meta_path = first_existing(meta_candidates)

    meta_summary = parse_salmon_meta(meta_path) if meta_path else {}
    total_reads = meta_summary.get("total_reads")
    mapped_reads = meta_summary.get("mapped_reads")

    if total_reads is not None:
        row["N reads/sample"] = total_reads
    if mapped_reads is not None:
        row["N uniquely mapped reads"] = mapped_reads

    uniq_frac = safe_fraction(mapped_reads, total_reads)
    if uniq_frac is not None:
        row["% uniquely mapped reads"] = convert_to_pct(uniq_frac)

    # ---- quant.json in *_counts → N cells, Total Genes Detected ----
    quant_candidates: List[Path] = []
    if counts_dir is not None:
        quant_candidates.extend(
            [
                counts_dir / "quant.json",
                counts_dir / "quant" / "quant.json",
            ]
        )
    quant_candidates.append(sample_root / "quant.json")
    quant_path = first_existing(quant_candidates)

    quant_summary = parse_alevin_quant(quant_path) if quant_path else {}

    if quant_summary.get("n_cells") is not None:
        row["N cells"] = quant_summary["n_cells"]
    if quant_summary.get("num_genes") is not None:
        row["Total Genes Detected"] = quant_summary["num_genes"]

    # ---- cell_meta.tsv → per-cell summaries (mean / medians) ----
    cell_meta_candidates: List[Path] = []
    if counts_dir is not None:
        cell_meta_candidates.extend(
            [
                counts_dir / "cell_meta.tsv",
                counts_dir / "cell_meta.txt",
                counts_dir / "quant" / "cell_meta.tsv",
            ]
        )
    if run_dir is not None:
        cell_meta_candidates.extend(
            [
                run_dir / "cell_meta.tsv",
                run_dir / "cell_meta.txt",
            ]
        )
    cell_meta_candidates.extend(
        [
            sample_root / "cell_meta.tsv",
            sample_root / "cell_meta.txt",
        ]
    )
    cell_meta_path = first_existing(cell_meta_candidates)

    if cell_meta_path:
        cell_summary = summarize_cell_meta(cell_meta_path)
    else:
        cell_summary = {}

    # Only fall back to cell_meta for N cells if quant.json didn't provide it
    if row["N cells"] is None and cell_summary.get("n_cells") is not None:
        row["N cells"] = cell_summary["n_cells"]

    if cell_summary.get("mean_reads_per_cell") is not None:
        row["Mean Reads per Cell"] = cell_summary["mean_reads_per_cell"]
    if cell_summary.get("median_umi_per_cell") is not None:
        row["Median UMI Counts per Cell"] = cell_summary["median_umi_per_cell"]
    if cell_summary.get("median_genes_per_cell") is not None:
        row["Median Genes per Cell"] = cell_summary["median_genes_per_cell"]

    return row


# ---------------------------------------------------------------------------
# Parse Biosciences helpers
# ---------------------------------------------------------------------------

def parse_parsebio_sample(sample_dir: Path) -> Dict[str, object]:
    """
    Parse one Parse Biosciences sample directory.
    """
    row: Dict[str, object] = {
        "Sample": sample_dir.name,
        "Software": "ParseBio_pipeline",
        "N reads/sample": None,
        "% uniquely mapped reads": None,
        "% multi-mapped reads": None,
        "N cells": None,
        "Mean Reads per Cell": None,
        "Median UMI Counts per Cell": None,
        "Median Genes per Cell": None,
    }

    report = sample_dir / "all-sample" / "report" / "sample_all_stats.csv"
    if not report.exists():
        print(f"[WARNING] ParseBio report not found: {report}")
        return row

    try:
        stats = pd.read_csv(
            report, header=None, names=["metric", "value"], index_col=0
        )["value"]

        total = clean_int(stats.get("reads_align_input"))
        unique = clean_int(stats.get("reads_align_unique"))
        multi = clean_int(stats.get("reads_align_multimap"))

        row["N reads/sample"] = total

        # Convert fractions to formatted percentages for ParseBio
        uniq_frac = safe_fraction(unique, total)
        multi_frac = safe_fraction(multi, total)

        row["% uniquely mapped reads"] = convert_to_pct(uniq_frac)
        row["% multi-mapped reads"] = convert_to_pct(multi_frac)

        row["N cells"] = clean_int(stats.get("number_of_cells"))

    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to parse ParseBio sample stats {report}: {e}")

    agg = sample_dir / "agg_samp_ana_summary.csv"
    if not agg.exists():
        print(f"[WARNING] ParseBio aggregate stats not found: {agg}")
        return row

    try:
        agg_stats = pd.read_csv(
            agg, header=None, names=["metric", "value"], index_col=0
        )["value"]
        row["Mean Reads per Cell"] = clean_int(agg_stats.get("mean_reads_per_cell"))
        row["Median UMI Counts per Cell"] = clean_int(
            agg_stats.get("ref-splitpipe_median_tscp_per_cell")
        )
        row["Median Genes per Cell"] = clean_int(
            agg_stats.get("ref-splitpipe_median_genes_per_cell")
        )
    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to parse ParseBio aggregate stats {agg}: {e}")

    return row


# ---------------------------------------------------------------------------
# Cell Ranger helpers
# ---------------------------------------------------------------------------

def parse_cellranger_sample(cr_dir: Path) -> Dict[str, object]:
    """
    Parse one Cell Ranger *_count result directory.
    """
    sample_name = cr_dir.name.replace("_count", "")
    row: Dict[str, object] = {
        "Sample": sample_name,
        "Software": "CellRanger_pipeline",
        "N cells": None,
        "Mean Reads per Cell": None,
        "Median Genes per Cell": None,
        "N reads/sample": None,
        "Saturation": None,
        "N R1 >Q30": None,
        "% uniquely mapped reads": None,
        "Total Genes Detected": None,
        "Median UMI Counts per Cell": None,
    }

    metrics_csv = cr_dir / "outs" / "log_metrics_summary.csv"
    if not metrics_csv.exists():
        print(f"[WARNING] Cell Ranger metrics CSV not found: {metrics_csv}")
        return row

    try:
        lines = metrics_csv.read_text().splitlines()
        if len(lines) < 2:
            print(f"[WARNING] Unexpected format in {metrics_csv} (too few lines)")
            return row

        second = re.sub(r"(\d),(\d)", r"\1\2", lines[1])
        cols = second.split(",")

        def col(idx: int, caster: Callable[[object], object] = clean_int):
            return caster(cols[idx]) if idx < len(cols) else None

        row["N cells"] = col(0)
        row["Mean Reads per Cell"] = col(1)
        row["Median Genes per Cell"] = col(2)
        row["N reads/sample"] = col(3)
        row["Saturation"] = col(6, clean_float)
        row["N R1 >Q30"] = col(8, clean_float)
        row["% uniquely mapped reads"] = col(11, clean_float)
        row["Total Genes Detected"] = col(18)
        row["Median UMI Counts per Cell"] = col(19)

    except Exception as e:  # noqa: BLE001
        print(f"[WARNING] Failed to parse Cell Ranger metrics {metrics_csv}: {e}")

    return row


# ---------------------------------------------------------------------------
# STARsolo traversal
# ---------------------------------------------------------------------------

def process_star_samples(
    star_root: Path,
    rows: List[Dict[str, object]],
    cellreads_dict: Dict[str, Dict[str, str]],
) -> None:
    """Fill rows and cellreads_dict with STARsolo samples."""
    for star_map_dir in sorted(star_root.glob("*")):
        if not star_map_dir.is_dir():
            continue

        print(f"Processing STARsolo sample: {star_map_dir.name}")
        sample = star_map_dir.name
        row: Dict[str, object] = {"Sample": sample, "Software": "STARsolo"}

        logfinal = find_star_file(star_map_dir, "Log.final.out")
        if logfinal:
            row.update(parse_star_logfinal(logfinal))

        solo_out = find_star_file(star_map_dir, "Solo.out")
        if solo_out is None:
            print(f"[WARNING] Solo.out not found in {star_map_dir}")
            rows.append(row)
            continue

        gene_dir, cfg_name = choose_gene_directory(solo_out)
        summary_df = load_summary_csv(gene_dir / "Summary.csv")
        row.update(extract_star_summary_fields(summary_df, cfg_name))

        filt_barcodes = gene_dir / "filtered" / "barcodes.tsv"
        row["N cells"] = count_file_lines(filt_barcodes)

        logout = find_star_file(star_map_dir, "Log.out")
        if logout:
            row["UMI cutoff used for cell calling"] = extract_umi_cutoff(
                logout, cfg_name
            )

        read_stats = gene_dir / "CellReads.stats"
        creads = add_cellreads_metrics(read_stats, filt_barcodes)
        cellreads_dict[sample] = creads

        creads_map: Mapping[str, str] = {
            "pct_noise": "Noise (% UMIs in non-cell barcodes)",
            "pct_exonic_reads": "% exonic reads",
            "pct_intronic_reads": "% Intronic reads",
            "pct_intergenic_reads": "% intergenic reads",
            "pct_mitochondrial_reads": "% mtDNA in Unique reads",
            "pct_exonicAS_reads": "% exonicAS reads",
            "pct_intronicAS_reads": "% intronicAS reads",
        }
        for key, col_name in creads_map.items():
            if key in creads:
                row[col_name] = creads[key]

        rows.append(row)


# ---------------------------------------------------------------------------
# TSV writer
# ---------------------------------------------------------------------------

TSV_FIELDS: List[str] = [
    "Sample",
    "Software",
    "N reads/sample",
    "Q30 Bases in CB+UMI",
    "Q30 Bases in RNA read",
    "N uniquely mapped reads",
    "% uniquely mapped reads",
    "% multi-mapped reads",
    "% multi-mapped reads: too many",
    "% unmapped: too short",
    "% unmapped: other",
    "Target N cells",
    "N cells",
    "UMI cutoff used for cell calling",
    "Saturation",
    "Reads for 0.7 saturation",
    "Noise (% UMIs in non-cell barcodes)",
    "% exonic reads",
    "% Intronic reads",
    "% intergenic reads",
    "% mtDNA in Unique reads",
    "% exonicAS reads",
    "% intronicAS reads",
    "% rRNA in Unique reads",
    "Mean Reads per Cell",
    "Median UMI Counts per Cell",
    "Median Genes per Cell",
    "Total Genes Detected",
]


def write_tsv_file(path: Path, rows_iter: Iterable[Mapping[str, object]]) -> None:
    """Write records to a TSV with canonical column order and blank NA cells."""
    df = pd.DataFrame(rows_iter)
    df = df.reindex(columns=TSV_FIELDS)
    df.to_csv(path, sep="\t", index=False, na_rep="")


# ---------------------------------------------------------------------------
# Main collection
# ---------------------------------------------------------------------------

def collect_all_log_metrics(
    root: Path,
) -> Tuple[List[Dict[str, object]], Dict[str, Dict[str, str]]]:
    """
    Collect per-sample QC rows and CellReads summaries from all supported
    pipelines under root.
    """
    summary_rows: List[Dict[str, object]] = []
    cellreads: Dict[str, Dict[str, str]] = {}

    # STARsolo
    star_root = root / "mapping_STARsolo"
    if star_root.exists():
        process_star_samples(star_root, summary_rows, cellreads)
    else:
        print(f"[WARNING] No STARsolo directory found at {star_root}")

    # Parse Biosciences
    parsebio_root = root / "ParseBio_pipeline"
    if parsebio_root.exists():
        for d in sorted(parsebio_root.glob("*")):
            if d.is_dir():
                summary_rows.append(parse_parsebio_sample(d))
    else:
        print(f"[WARNING] No ParseBio directory found at {parsebio_root}")

    # Cell Ranger
    cellranger_root = root / "CellRanger_pipeline"
    if cellranger_root.exists():
        for d in sorted(cellranger_root.glob("*")):
            if d.is_dir():
                summary_rows.append(parse_cellranger_sample(d))
    else:
        print(f"[WARNING] No CellRanger directory found at {cellranger_root}")

    # sci-rocket
    sci_root = root / "sci-rocket"
    if sci_root.exists():
        qc_files = list(sci_root.glob("*_scirocket/sci-dash/js/qc_data.js"))
        if not qc_files:
            print(f"[WARNING] No sci-rocket qc_data.js found under {sci_root}")
        for qc in sorted(qc_files):
            print(f"Processing sci-rocket QC: {qc}")
            summary_rows.extend(parse_scirocket_run(qc))
    else:
        print(f"[WARNING] No sci-rocket directory found at {sci_root}")

    # alevin-fry + alevinQC (mapping_alevin)
    alevin_root = root / "mapping_alevin"
    if alevin_root.exists():
        for sample_root in sorted(alevin_root.glob("*")):
            if sample_root.is_dir():
                print(f"Processing alevin-fry sample: {sample_root.name}")
                summary_rows.append(parse_alevinfry_sample(sample_root))
    else:
        print(f"[WARNING] No alevin-fry directory found at {alevin_root}")

    return summary_rows, cellreads


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> None:
    """CLI entry point. Writes TSV and optional JSON summary."""
    args = parse_command_line_arguments(argv)

    rows, cellreads = collect_all_log_metrics(args.outdir)
    write_tsv_file(args.outfile, rows)

    if not rows:
        print(f"[WARNING] No mapping stats found in {args.outdir}")
        write_tsv_file(args.outfile, [])

    if args.json:
        try:
            args.json.write_text(json.dumps(cellreads, indent=2))
            print(f"Wrote CellReads JSON → {args.json}")
        except Exception as e:  # noqa: BLE001
            print(f"[WARNING] Could not write JSON to {args.json}: {e}")

    print(f"Done. TSV rows written: {len(rows)}")


if __name__ == "__main__":
    main()
