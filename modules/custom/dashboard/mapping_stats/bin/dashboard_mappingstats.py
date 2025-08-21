#!/usr/bin/env python
"""

Aggregate QC log metrics from STARsolo, ParseBio and CellRanger pipelines
into a single TSV with universal column names. Optionally writes a JSON
file with CellReads-derived percentages per sample.

Example
-------
Run on a results directory that contains subfolders like:
    mapping_STARsolo/, ParseBio_pipeline/, CellRanger_pipeline/

    python dashboard_mappingstats.py /analysis/ \\
        --outfile mapping_stats.tsv \\
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

def parse_command_line_arguments(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    Returns:
        Parsed arguments namespace with:
          - output_dir (Path): Root directory with pipeline results.
          - outfile (Path): Destination TSV path.
          - json (Optional[Path]): Optional JSON output path for CellReads summaries.
    """
    parser = argparse.ArgumentParser(
        description="Aggregate mapping QC log_metrics into one TSV file."
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Results directory containing mapping_STARsolo/, ParseBio_pipeline/, CellRanger_pipeline/ …",
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
        help="Optional JSON with CellReads-derived percentages",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def run_command(cmd: List[str]) -> str:
    """
    Run a shell command and return stdout.
    On failure (non-zero exit or missing binary) prints a warning and returns an empty string.
    """
    try:
        return subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"[WARNING] Command failed: {' '.join(cmd)} -> {e}")
        return ""


def clean_int(val: object) -> Optional[int]:
    """
    Parse an integer from a string that may contain commas or a trailing '%'.
    """
    try:
        return int(str(val).replace(",", "").replace("%", ""))
    except (TypeError, ValueError):
        return None


def clean_float(val: object) -> Optional[float]:
    """
    Parse a float from a string that may contain a trailing '%'.
    """
    try:
        return float(str(val).replace("%", ""))
    except (TypeError, ValueError):
        return None


def convert_to_pct(value: Optional[float]) -> str:
    """
    Format a fraction (0..1) as a percentage string with two decimals.
    """
    if value is None:
        return ""
    return f"{value * 100:.2f}%"


# ---------------------------------------------------------------------------
# STARsolo helpers
# ---------------------------------------------------------------------------

# Map of metric names
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
    """
    Extract selected metrics from STARsolo `Log.final.out`.
    """
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
    except Exception as e:
        print(f"[WARNING] Failed to parse {log_file}: {e}")
    return metrics


def choose_gene_directory(solo_out: Path) -> Tuple[Path, str]:
    """
    Select the STARsolo gene results directory and its config label.
    Prefers `GeneFull_Ex50pAS` when present; otherwise returns `Gene`.
    """
    gene_full = solo_out / "GeneFull_Ex50pAS"
    return (gene_full, "GeneFull_Ex50pAS") if gene_full.exists() else (solo_out / "Gene", "Gene")


def load_summary_csv(csv_path: Path) -> pd.DataFrame:
    """
    Load STARsolo Summary.csv as a 2-column table with the first column as index.
    """
    if not csv_path.exists():
        print(f"[WARNING] Summary CSV file not found: {csv_path}")
        return pd.DataFrame()
    try:
        df = pd.read_csv(csv_path)
        df.set_index(df.columns[0], inplace=True)
        return df
    except Exception as e:
        print(f"[WARNING] Could not read STARsolo Summary.csv at {csv_path}: {e}")
        return pd.DataFrame()


def extract_star_summary_fields(df: pd.DataFrame, cfg_name: str) -> Dict[str, str]:
    """
    Pick and normalize selected fields from STARsolo Summary.csv.
    Metrics that represent proportions (Q30, Saturation) are formatted as percentages.
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
        if value in df.index:
            try:
                val = df.at[value, df.columns[-1]]
                out[key] = convert_to_pct(val) if ("Q30" in key or "Saturation" in key) else str(int(val))
            except Exception:
                out[key] = ""
    return out


def count_file_lines(p: Path) -> Optional[int]:
    """
    Count number of lines in a text file.
    """
    if not p.exists():
        print(f"[WARNING] File not found (cannot count lines): {p}")
        return None
    try:
        with p.open("rb") as fh:
            return sum(1 for _ in fh)
    except Exception as e:
        print(f"[WARNING] Failed to count lines in {p}: {e}")
        return None


def extract_umi_cutoff(log_out: Path, cfg_name: str) -> Optional[int]:
    """
    Extract `nUMImin` (cell-calling cutoff) from STARsolo Log.out.
    Scans the section corresponding to the given cfg (Gene/GeneFull_Ex50pAS).
    """
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
                m = re.search(r"nUMImin=(\d+)", line)
                if m:
                    return int(m.group(1))
    except Exception as e:
        print(f"[WARNING] Failed to parse {log_out}: {e}")
    return None


# ---------------------------------------------------------------------------
# CellReads helpers  (STARsolo)
# ---------------------------------------------------------------------------

def add_cellreads_metrics(read_stats: Path, barcodes: Path) -> Dict[str, str]:
    """
    Compute CellReads-derived percentages as strings with '%' suffix.
    Returns dict with keys like 'pct_noise', 'pct_exonic_reads', ..., each as '12.34%'.
    """
    if not read_stats.exists() or not barcodes.exists():
        print(f"[WARNING] CellReads stats or barcodes file not found: {read_stats}, {barcodes}")
        return {}

    try:
        stats_df = pd.read_csv(read_stats, sep="\t", index_col=0)
        filtered_cbs = set(pd.read_csv(barcodes, header=None)[0])

        total_cbmatch = stats_df.get("cbMatch", pd.Series(dtype=float)).sum()
        cbmatch_noise = stats_df.loc["CBnotInPasslist", "cbMatch"] if "CBnotInPasslist" in stats_df.index else 0

        # Restrict to passing barcodes for genomic composition
        filtered_df = stats_df.loc[stats_df.index.intersection(filtered_cbs)]
        summed_fil = filtered_df.sum(numeric_only=True)

        genome_total = summed_fil.get("genomeU", 0) + summed_fil.get("genomeM", 0) or None

        out: Dict[str, str] = {}

        # Noise (% reads/UMIs with CB not in pass-list)
        if total_cbmatch:
            out["pct_noise"] = convert_to_pct(cbmatch_noise / total_cbmatch)

        # Genomic composition (% of reads mapped to categories)
        if genome_total:
            out["pct_exonic_reads"]        = convert_to_pct(summed_fil.get("exonic",    0) / genome_total)
            out["pct_intronic_reads"]      = convert_to_pct(summed_fil.get("intronic",  0) / genome_total)
            out["pct_intergenic_reads"]    = convert_to_pct(
                (genome_total
                 - summed_fil.get("exonic",    0)
                 - summed_fil.get("intronic",  0)
                 - summed_fil.get("exonicAS",  0)
                 - summed_fil.get("intronicAS",0)
                ) / genome_total
            )
            out["pct_mitochondrial_reads"] = convert_to_pct(summed_fil.get("mito",      0) / genome_total)
            out["pct_exonicAS_reads"]      = convert_to_pct(summed_fil.get("exonicAS",  0) / genome_total)
            out["pct_intronicAS_reads"]    = convert_to_pct(summed_fil.get("intronicAS",0) / genome_total)

        print("Extracted CellReads metrics:", out)
        return out
    
    except Exception as e:
        print(f"[WARNING] Failed to compute CellReads metrics from {read_stats}: {e}")
        return {}


# ---------------------------------------------------------------------------
# ParseBio helpers
# ---------------------------------------------------------------------------

def parse_parsebio_sample(sample_dir: Path) -> Dict[str, object]:
    """
    Parses one ParseBio sample directory.
    Returns a dict keyed by TSV column names for each sample.
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
        s = pd.read_csv(report, header=None, names=["metric", "value"], index_col=0)["value"]
        total = clean_int(s.get("reads_align_input"))
        unique = clean_int(s.get("reads_align_unique"))
        multi = clean_int(s.get("reads_align_multimap"))

        row.update({
            "N reads/sample": total,
            "% uniquely mapped reads": (unique / total) if total else None,
            "% multi-mapped reads": (multi / total) if total else None,
            "N cells": clean_int(s.get("number_of_cells")),
        })
    except Exception as e:
        print(f"[WARNING] Failed to parse ParseBio sample stats {report}: {e}")

    agg = sample_dir / "agg_samp_ana_summary.csv"
    if agg.exists():
        try:
            agg_s = pd.read_csv(agg, header=None, names=["metric", "value"], index_col=0)["value"]
            row.update({
                "Mean Reads per Cell": clean_int(agg_s.get("mean_reads_per_cell")),
                "Median UMI Counts per Cell": clean_int(agg_s.get("ref-splitpipe_median_tscp_per_cell")),
                "Median Genes per Cell": clean_int(agg_s.get("ref-splitpipe_median_genes_per_cell")),
            })
        except Exception as e:
            print(f"[WARNING] Failed to parse ParseBio aggregate stats {agg}: {e}")
    else:
        print(f"[WARNING] ParseBio aggregate stats not found: {agg}")

    return row


# ---------------------------------------------------------------------------
# Cell Ranger helpers
# ---------------------------------------------------------------------------

def parse_cellranger_sample(cr_dir: Path) -> Dict[str, object]:
    """
    Parses one Cell Ranger `*_count` result directory.
    Returns a dict keyed by TSV column names for each sample.
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
        # Commas inside thousands need stripping on the second line
        second = re.sub(r"(\d),(\d)", r"\1\2", lines[1])
        cols = second.split(",")

        def col(i: int, caster: Callable[[object], object] = clean_int):
            return caster(cols[i]) if i < len(cols) else None

        row.update({
            "N cells": col(0),
            "Mean Reads per Cell": col(1),
            "Median Genes per Cell": col(2),
            "N reads/sample": col(3),
            "Saturation": col(6, clean_float),
            "N R1 >Q30": col(8, clean_float),
            "% uniquely mapped reads": col(11, clean_float),
            "Total Genes Detected": col(18),
            "Median UMI Counts per Cell": col(19),
        })
    except Exception as e:
        print(f"[WARNING] Failed to parse Cell Ranger metrics {metrics_csv}: {e}")

    return row


# ---------------------------------------------------------------------------
# Traversal / aggregation
# ---------------------------------------------------------------------------

def find_star_file(sample_dir: Path, suffix: str) -> Optional[Path]:
    """
    Find a STARsolo file (Log.final.out, Log.out, Solo.out directory) that may have
    a sample prefix (e.g. r12345_Log.final.out).
    Returns the first match or None.
    """
    # Exact match first
    p = sample_dir / suffix
    if p.exists():
        return p

    # Glob for *_suffix
    matches = list(sample_dir.glob(f"*_{suffix}"))
    if matches:
        return matches[0]

    return None



def process_star_samples(star_root: Path,
                         rows: List[Dict[str, object]],
                         cellreads_dict: Dict[str, Dict[str, str]]) -> None:
    for star_map_dir in sorted(star_root.glob("*")):
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
            row["UMI cutoff used for cell calling"] = extract_umi_cutoff(logout, cfg_name)

        # CellReads percentages
        read_stats = gene_dir / "CellReads.stats"
        creads = add_cellreads_metrics(read_stats, filt_barcodes)
        cellreads_dict[sample] = creads  # store even if empty

        # Map CellReads keys to TSV column names (only present keys are added)
        CREADS_MAP: Mapping[str, str] = {
            "pct_noise": "Noise (% UMIs in non-cell barcodes)",
            "pct_exonic_reads": "% exonic reads",
            "pct_intronic_reads": "% Intronic reads",
            "pct_intergenic_reads": "% intergenic reads",
            "pct_mitochondrial_reads": "% mtDNA in Unique reads",
            "pct_exonicAS_reads": "% exonicAS reads",
            "pct_intronicAS_reads": "% intronicAS reads",
        }
        for metric_key, col_name in CREADS_MAP.items():
            if metric_key in creads:
                row[col_name] = creads[metric_key]

        rows.append(row)  # always append, even if mostly empty


# ---------------------------------------------------------------------------
# TSV writer
# ---------------------------------------------------------------------------

#: Canonical TSV column order.
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
    """
    Write records to a TSV with canonical column order and blank NA cells.
    """
    df = pd.DataFrame(rows_iter)
    # Ensure every expected column is present and ordered
    df = df.reindex(columns=TSV_FIELDS)
    df.to_csv(path, sep="\t", index=False, na_rep="")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def collect_all_log_metrics(root: Path) -> Tuple[List[Dict[str, object]], Dict[str, Dict[str, str]]]:
    """
    Searches root directory and collect per-sample QC rows and CellReads summaries.
    """
    summary_rows: List[Dict[str, object]] = []
    cellreads: Dict[str, Dict[str, str]] = {}

    # STARsolo
    star_root = root / "mapping_STARsolo"
    if star_root.exists():
        process_star_samples(star_root, summary_rows, cellreads)
    else:
        print(f"[WARNING] No STARsolo directory found at {star_root}")

    # ParseBio
    parsebio_root = root / "ParseBio_pipeline"
    if parsebio_root.exists():
        for d in sorted(parsebio_root.glob("*")):
            summary_rows.append(parse_parsebio_sample(d))
    else:
        print(f"[WARNING] No ParseBio directory found at {parsebio_root}")

    # Cell Ranger
    cellranger_root = root / "CellRanger_pipeline"
    if cellranger_root.exists():
        for d in sorted(cellranger_root.glob("*")):
            summary_rows.append(parse_cellranger_sample(d))
    else:
        print(f"[WARNING] No CellRanger directory found at {cellranger_root}")

    return summary_rows, cellreads


def main(argv: Optional[List[str]] = None) -> None:
    """
    Writes the TSV and optional JSON.
    Prints warnings for missing inputs; never raises error due to missing files.
    """
    args = parse_command_line_arguments(argv)

    rows, cellreads = collect_all_log_metrics(args.output_dir)
    write_tsv_file(args.outfile, rows)

    if not rows:
        print(f"[WARNING] No mapping stats found in {args.output_dir}")
        # still ensure an empty TSV with headers exists
        write_tsv_file(args.outfile, [])

    if args.json:
        try:
            args.json.write_text(json.dumps(cellreads, indent=2))
            print(f"Wrote CellReads JSON → {args.json}")
        except Exception as e:
            print(f"[WARNING] Could not write JSON to {args.json}: {e}")

    print(f"Done. TSV rows written: {len(rows)}")


if __name__ == "__main__":
    main()
