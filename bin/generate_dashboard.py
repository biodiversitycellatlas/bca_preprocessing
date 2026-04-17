#!/usr/bin/env python3
"""
Embed data into HTML dashboard summarizing scRNA-seq preprocessing run.
"""

import argparse
import json
import csv
import base64
import re
import os
import sys
import datetime
from typing import Dict, Any, Optional, List

def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description="Generate HTML Dashboard")
    parser.add_argument("--template", required=True, help="Path to HTML template")
    parser.add_argument("--output", required=True, help="Output HTML filename")

    # Core Metadata
    parser.add_argument("--project", default="Biodiversity Cell Atlas", help="Project Name")
    parser.add_argument("--pipeline", default="bca-preprocessing", help="Pipeline Name")
    parser.add_argument("--version", default="0.2.1", help="Pipeline Version")
    parser.add_argument("--commit", default="unknown", help="Git Commit Hash")

    # Configurations & Manifest
    parser.add_argument("--run_config", help="Path to run_config_date.txt")
    parser.add_argument("--samplesheet", required=True, help="Path to samplesheet.csv")
    parser.add_argument("--analytical_samples", required=True, help="CSV manifest from Nextflow")

    # File lists
    # Data Sources (STAR / STARsolo)
    parser.add_argument("--star_logs", nargs='*', help="List of STAR Log.final.out files")
    parser.add_argument("--star_summaries", nargs='*', help="List of STARsolo Summary.csv files")
    parser.add_argument("--star_full_logs", nargs='*', help="List of STAR Log.out files (for UMI threshold)")
    parser.add_argument("--mt_rrna_metrics", nargs='*', help="List of *_mt_rrna_metrics.txt files")
    parser.add_argument("--saturation_logs", nargs='*', help="List of saturation.log files")
    parser.add_argument("--cell_stats", nargs='*', help="List of STARsolo CellReads.stats files")

    # Data Sources (alevin-fry fallback; optional)
    parser.add_argument("--af_meta_info", nargs='*', help="List of meta_info.json files (e.g. total_reads, mapping_rate)")
    parser.add_argument("--af_quant_json", nargs='*', help="List of quant.json files from `alevin-fry quant` output")
    parser.add_argument("--af_cell_meta", nargs='*', help="List of cell_meta.tsv files (per-cell metrics; mean/median computed when possible)")

    # Visualization Files
    parser.add_argument("--sankey_files", nargs='*', help="List of *_kraken.sankey.html files")
    parser.add_argument("--per_cell_files", nargs='*', help="List of per-cell JSON files (e.g. *_metrics.json)")
    parser.add_argument("--saturation_imgs", nargs='*', help="List of saturation PNG images")
    parser.add_argument("--residuals_imgs", nargs='*', help="List of residuals PNG images")
    parser.add_argument("--knee_files", nargs='*', help="List of UMIperCellSorted.txt files")

    parser.add_argument("--cellsweep_tables", nargs='*')
    parser.add_argument("--cellsweep_plots_contrib", nargs='*')
    parser.add_argument("--cellsweep_plots_umap", nargs='*')
    parser.add_argument("--cellbender_tables", nargs='*')
    parser.add_argument("--cellbender_plots1", nargs='*')
    parser.add_argument("--cellbender_plots2", nargs='*')

    return parser.parse_args()


# --- Formatters & Encoders ---

def to_pct(val: Any, precision: int = 2) -> str:
    """
    Convert a fraction (0..1) to a percentage string.
    Returns "N/A" for None/empty/"N/A"/unparseable values.
    """
    if val is None or str(val).strip().upper() == "N/A" or val == "":
        return "N/A"
    try:
        f = float(val)
        return f"{f * 100:.{precision}f}%"
    except (ValueError, TypeError):
        return "N/A"


def encode_image(image_path: Optional[str]) -> Optional[str]:
    """
    Return a PNG data URI for embedding in HTML, or None if missing/unreadable.
    """
    if not image_path or not os.path.exists(image_path):
        return None
    try:
        with open(image_path, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
        return f"data:image/png;base64,{encoded_string}"
    except Exception as e:
        sys.stderr.write(f"Error encoding image {image_path}: {e}\n")
        return None

def safe_float(val: Any) -> Optional[float]:
    try: return float(val)
    except (ValueError, TypeError): return None

def fmt(val: Any) -> str:
    if val is None or val == "": return ""
    try:
        s = str(val).strip()
        if s.isdigit(): return "{:,}".format(int(s))
    except Exception: pass
    return str(val)


# --- Parsers ---

def extract_sankey_data(html_path: Optional[str]) -> Optional[dict]:
    """
    Extract htmlwidgets JSON from a sankey HTML file, or None if missing/unparseable.
    """
    if not html_path or not os.path.exists(html_path):
        return None
    try:
        with open(html_path, 'r') as f: content = f.read()
        match = re.search(r'<script type="application/json" data-for="htmlwidget-[a-f0-9]+">(.*?)</script>', content, re.DOTALL)
        if match: return json.loads(match.group(1)).get('x', None)
    except Exception: pass
    return None

def parse_saturation_log(path: Optional[str]) -> str:
    """
    Extract read input needed to achieve n saturation from saturation.log.
    Returns a string like '22.5 M' or 'N/A'.
    """
    if not path or not os.path.exists(path):
        return "N/A"
    try:
        with open(path, 'r') as f: content = f.read()
        match = re.search(r"To achieve a saturation of [0-9.]+, ninput should be approximately:\s*(.*?)\s*reads", content)
        if match: return match.group(1).strip()
    except Exception: pass
    return "N/A"

def parse_run_config(path: Optional[str]) -> Dict[str, str]:
    """
    Parse a KEY=VALUE text file into a dict; return {} if missing/unreadable.
    """
    config: Dict[str, str] = {}
    if not path or not os.path.exists(path):
        return config
    try:
        with open(path, 'r') as f:
            for line in f:
                if '=' in line:
                    k, v = line.split('=', 1)
                    config[k.strip()] = v.strip()
    except Exception: pass
    return config

def parse_cell_filtering_table(path: Optional[str], method: str) -> List[Dict[str, Any]]:
    if not path or not os.path.exists(path): return []
    data = []
    try:
        with open(path, 'r', newline='', encoding='utf-8-sig') as f:  # utf-8-sig strips BOM
            header_line = f.readline()
            delimiter = '\t' if '\t' in header_line else ','
            f.seek(0)
            reader = csv.DictReader(f, delimiter=delimiter)
            for row in reader:
                cleaned_row = {k.strip(): v.strip() for k, v in row.items() if k}  # drop `and v`
                if method == "cellsweep":
                    data.append({
                        "gene": cleaned_row.get("gene", "N/A"),
                        "ambient_hat": safe_float(cleaned_row.get("ambient_hat"))
                    })
                elif method == "cellbender":
                    data.append({
                        "gene_name": cleaned_row.get("gene_name", "N/A"),
                        "n_removed": safe_float(cleaned_row.get("n_removed")),
                        "fraction_removed": safe_float(cleaned_row.get("fraction_removed"))
                    })
        if method == "cellsweep":
            data.sort(key=lambda x: x["ambient_hat"] if x["ambient_hat"] is not None else -1, reverse=True)
        elif method == "cellbender":
            data.sort(key=lambda x: x["n_removed"] if x["n_removed"] is not None else -1, reverse=True)
        return data[:100]
    except Exception: return []

def parse_star_log(path: Optional[str]) -> Dict[str, str]:
    """Parse STAR Log.final.out (key | value) into a dict."""
    data: Dict[str, str] = {}
    if not path or not os.path.exists(path):
        return data
    try:
        with open(path, 'r') as f:
            for line in f:
                if '|' in line:
                    parts = line.split('|', 1)
                    data[parts[0].strip()] = parts[1].strip()
    except Exception: pass
    return data

def parse_starsolo_summary(path: Optional[str]) -> Dict[str, str]:
    """
    Parse STARsolo Summary.csv into a dict.
    """
    data: Dict[str, str] = {}
    if not path or not os.path.exists(path):
        return data
    try:
        with open(path, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) >= 2: data[row[0].strip()] = row[1].strip()
    except Exception: pass
    return data

def parse_mt_rrna_metrics(path: Optional[str]) -> Dict[str, str]:
    """
    Parse MT/rRNA metrics file with 'key,value' lines into a dict.
    """
    data: Dict[str, str] = {}
    if not path or not os.path.exists(path):
        return data
    try:
        with open(path, 'r') as f:
            for line in f:
                if ',' in line:
                    k, v = line.strip().split(',', 1)
                    data[k.strip()] = v.strip()
    except Exception: pass
    return data

def parse_starsolo_intronic(path: Optional[str]) -> Any:
    """
    Compute intronic fraction from STARsolo CellReads.stats.
    Returns: float fraction or "N/A"
    """
    if not path or not os.path.exists(path):
        return "N/A"

    try:
        mapped_reads = 0
        intronic_sum = 0
        with open(path, 'r') as f:
            next(f)
            for line in f:
                parts = line.strip().split()
                if len(parts) < 6: continue
                intronic_sum += (int(parts[4]) + int(parts[5]))
                mapped_reads += (int(parts[2]) + int(parts[3]) + int(parts[4]) + int(parts[5]))
        if mapped_reads > 0: return intronic_sum / mapped_reads
    except Exception: pass
    return "N/A"

def extract_umi_cutoff(log_out_path: Optional[str], cfg_name: str = "GeneFull_Ex50pAS") -> int:
    """
    Extract nUMImin from STAR Log.out within the Solo post-map block.
    Returns 0 if unavailable.
    """
    if not log_out_path or not os.path.exists(log_out_path):
        return 0
    try:
        with open(log_out_path, 'r') as f:
            in_block = False
            for line in f:
                if f"Starting Solo post-map for {cfg_name}" in line: in_block = True
                elif in_block and "Starting Solo post-map for" in line: break
                elif in_block and "cellFiltering" in line:
                    match = re.search(r"nUMImin=(\d+)", line)
                    if match: return int(match.group(1))
    except Exception: pass
    return 0

def _safe_read_json(path: Optional[str]) -> Dict[str, Any]:
    """
    Read a JSON object; return {} if missing/unreadable/non-dict.
    """
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path, "r") as f:
            obj = json.load(f)
        return obj if isinstance(obj, dict) else {}
    except Exception: return {}

def parse_cell_meta_tsv(path: Optional[str]) -> Dict[str, Any]:
    """
    Parse a per-cell TSV and compute basic per-column mean/median for numeric columns.
    """
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path, "r", newline="") as f:
            rows = [r for r in csv.DictReader(f, delimiter="\t") if r]
        if not rows: return {}

        numeric = {}
        for r in rows:
            for k, v in r.items():
                if v is None or str(v).strip() == "": continue
                try: numeric.setdefault(k, []).append(float(v))
                except ValueError: pass

        def median(xs):
            ys = sorted(xs)
            mid = len(ys) // 2
            return ys[mid] if len(ys) % 2 == 1 else (ys[mid - 1] + ys[mid]) / 2.0

        out = {"n_rows": len(rows)}
        for k, xs in numeric.items():
            if len(xs) >= 2:
                out[f"{k}__mean"] = sum(xs) / len(xs)
                out[f"{k}__median"] = median(xs)
        return out
    except Exception: return {}


def main() -> None:
    args = parse_args()

    # 1) Load run config
    raw_run_config = parse_run_config(args.run_config)

    # 2) Report metadata
    report_metadata = {
        "project": args.project,
        "pipeline": args.pipeline,
        "version": args.version,
        "date": datetime.date.today().isoformat(),
        "commit": args.commit
    }

    samplesheet_config = {}
    if os.path.exists(args.samplesheet):
        with open(args.samplesheet, 'r') as f:
            for row in csv.DictReader(f):
                if row.get('sample'): samplesheet_config[row['sample']] = row

    # Load Nextflow Manifest mapping analytical IDs to Base IDs
    active_samples = []
    if os.path.exists(args.analytical_samples):
        with open(args.analytical_samples, 'r') as f:
            for row in csv.DictReader(f):
                active_samples.append({
                    "id": row["analytical_id"],
                    "base": row["base_id"],
                    "source": row.get("source", "auto").strip().lower()
                })
    else:
        sys.exit(f"Error: Manifest {args.analytical_samples} not found.")

    print(f"Loaded {len(active_samples)} analytical samples from Nextflow manifest.", file=sys.stderr)

    # Localized fast file mapping, sort IDs descending by length so 'sample1_subsampled_starsolo'
    file_map = {s["id"]: {} for s in active_samples}
    sample_info_by_id = {s["id"]: s for s in active_samples}
    sorted_ids = sorted([s["id"] for s in active_samples], key=len, reverse=True)

    def map_files(cli_list: Optional[List[str]], key_name: str, allowed_sources: Optional[set] = None):
        if not cli_list:
            return
        for f in cli_list:
            if not f:
                continue
            base_name = os.path.basename(f)
            for sid in sorted_ids:
                if sid not in base_name:
                    continue
                src = sample_info_by_id[sid]["source"]
                if allowed_sources and src not in allowed_sources:
                    continue
                file_map[sid][key_name] = f
                break

    # Keep shared only if both pipelines can emit it, otherwise assign to the one that emits it
    map_files(args.star_logs, 'star_log', {'starsolo'})
    map_files(args.star_summaries, 'star_summary', {'starsolo'})
    map_files(args.star_full_logs, 'star_full_log', {'starsolo'})
    map_files(args.mt_rrna_metrics, 'mt_rrna', {'starsolo'})
    map_files(args.saturation_logs, 'sat_log', {'starsolo'})
    map_files(args.cell_stats, 'cell_stats', {'starsolo'})
    map_files(args.saturation_imgs, 'sat_img', {'starsolo'})
    map_files(args.residuals_imgs, 'res_img', {'starsolo'})
    map_files(args.sankey_files, 'sankey', {'starsolo'})
    map_files(args.per_cell_files, 'per_cell', {'starsolo'})
    map_files(args.knee_files, 'knee', {'starsolo'})

    map_files(args.af_meta_info, 'af_meta', {'alevin'})
    map_files(args.af_quant_json, 'af_quant', {'alevin'})
    map_files(args.af_cell_meta, 'af_cell', {'alevin'})

    # Shared exception
    map_files(args.cellsweep_tables, 'cs_table')
    map_files(args.cellsweep_plots_contrib, 'cs_plot_contrib')
    map_files(args.cellsweep_plots_umap, 'cs_plot_umap')
    map_files(args.cellbender_tables, 'cb_table')
    map_files(args.cellbender_plots1, 'cb_plot1')
    map_files(args.cellbender_plots2, 'cb_plot2')


    # Output structures
    global_cols = [
        "Sample", "% Uniquely Mapped Reads", "N cells", "Saturation",
        "Reads Needed for Target Saturation", "Noise (% UMIs non-cell barcodes)",
        "Median Transcripts Per Cell", "% Intronic Reads", "% rRNA in Unique reads",
        "% mtDNA in Unique reads", "% mtDNA in multimappers all pos"
    ]
    global_rows = []
    samples_json_list = []
    per_cell_data = {}
    saturation_images = {}
    knee_data = {}
    cell_filtering_data = {}

    for sample_info in active_samples:
        s_id = sample_info["id"]
        base_id = sample_info["base"]
        files = file_map[s_id]

        print(f"Processing {s_id}...", file=sys.stderr)

        source = sample_info.get("source", "auto")
        using_star = source == "starsolo" if source != "auto" else bool(files.get('star_log') and files.get('star_summary'))
        using_alevin = source == "alevin" if source != "auto" else not using_star

        mt_stats = parse_mt_rrna_metrics(files.get('mt_rrna'))
        reads_07_sat_val = parse_saturation_log(files.get('sat_log')) if using_star else "N/A"
        intronic_pct = to_pct(parse_starsolo_intronic(files.get('cell_stats'))) if using_star else "N/A"

        umi_threshold = 0
        if using_star:
            full_log = files.get('star_full_log')
            if not full_log and files.get('star_log'):
                guess_1 = files['star_log'].replace("Log.final.out", "Log.out")
                if os.path.exists(guess_1): full_log = guess_1
            umi_threshold = extract_umi_cutoff(full_log)

        def get_val(source, key, default="N/A"):
            return source.get(key, default) if source else default

        # Metrics Assembly
        if using_star:
            star_stats = parse_star_log(files.get('star_log'))
            solo_stats = parse_starsolo_summary(files.get('star_summary'))

            n_input_reads = get_val(star_stats, "Number of input reads")
            n_unique = get_val(star_stats, "Uniquely mapped reads number")
            pct_unique = get_val(star_stats, "Uniquely mapped reads %")
            pct_multi = get_val(star_stats, "Number of reads mapped to multiple loci")
            if "%" not in str(pct_multi) and n_input_reads != "N/A" and pct_multi != "N/A":
                pct_multi = get_val(star_stats, "% of reads mapped to multiple loci", pct_multi)

            pct_multi_too_many = get_val(star_stats, "% of reads mapped to too many loci")
            pct_short = get_val(star_stats, "% of reads unmapped: too short")
            pct_other = get_val(star_stats, "% of reads unmapped: other")
            pct_r1_q30 = to_pct(get_val(solo_stats, "Q30 Bases in CB+UMI"))
            pct_r2_q30 = to_pct(get_val(solo_stats, "Q30 Bases in RNA read"))
            saturation = to_pct(get_val(solo_stats, "Sequencing Saturation"))

            n_cells = get_val(solo_stats, "Estimated Number of Cells", get_val(solo_stats, "Cells Detected"))
            mean_reads = get_val(solo_stats, "Mean Reads per Cell")
            median_umis = get_val(solo_stats, "Median UMI per Cell")
            median_genes = get_val(solo_stats, "Median GeneFull_Ex50pAS per Cell")
            total_genes = get_val(solo_stats, "Total GeneFull_Ex50pAS Detected")

            frac_in_cells = get_val(solo_stats, "Fraction of Unique Reads in Cells")
            noise_pct = "N/A"
            if frac_in_cells != "N/A":
                try: noise_pct = f"{(1.0 - float(frac_in_cells)) * 100:.2f}%"
                except: pass

        else:
            af_meta = _safe_read_json(files.get('af_meta'))
            af_quant = _safe_read_json(files.get('af_quant'))
            af_cell = parse_cell_meta_tsv(files.get('af_cell'))

            total_reads = af_meta.get("total_reads", af_meta.get("num_processed", "N/A"))
            mapping_rate = af_meta.get("mapping_rate", "N/A")
            if mapping_rate == "N/A" and "percent_mapped" in af_meta:
                mapping_rate = float(af_meta["percent_mapped"]) / 100.0

            mapped_reads = af_meta.get("num_mapped", "N/A")
            if mapped_reads == "N/A" and total_reads != "N/A" and mapping_rate != "N/A":
                try: mapped_reads = int(float(total_reads) * float(mapping_rate))
                except: pass

            n_input_reads = total_reads
            n_unique = mapped_reads
            pct_unique = to_pct(mapping_rate)
            n_cells = af_quant.get("num_quantified_cells", af_quant.get("num_cells", af_meta.get("final_num_cbs", "N/A")))
            total_genes = af_quant.get("num_genes", af_quant.get("num_targets", "N/A"))

            mean_reads = "N/A"
            if total_reads != "N/A" and n_cells != "N/A" and int(n_cells) > 0:
                try: mean_reads = int(float(total_reads) / int(n_cells))
                except: pass

            pct_multi = pct_multi_too_many = pct_short = pct_other = "N/A"
            pct_r1_q30 = pct_r2_q30 = saturation = noise_pct = "N/A"

            median_umis = af_cell.get("nUMI__median", af_cell.get("umis__median", af_meta.get("mean_umis_per_cell", "N/A")))
            median_genes = af_cell.get("nGene__median", af_cell.get("genes__median", af_meta.get("mean_genes_per_cell", "N/A")))

        rrna_pct = to_pct(get_val(mt_stats, "Percentage of rRNA reads (of uniquely mapped reads)"))
        mtdna_unique = to_pct(get_val(mt_stats, "Percentage of mtDNA reads (of mapped reads)"))
        mtdna_multi_all = to_pct(get_val(mt_stats, "Percentage of mtDNA in multimapped reads (all alignments)"))
        mtdna_multi_primary = to_pct(get_val(mt_stats, "Percentage of mtDNA in multimapped reads (primary alignment)"))

        global_rows.append([
            s_id, pct_unique, fmt(n_cells), saturation, reads_07_sat_val, noise_pct,
            fmt(median_umis), intronic_pct, rrna_pct, mtdna_unique, mtdna_multi_all
        ])

        # Cell Filtering Data
        cf_dict = None
        if files.get('cs_table'):
            cf_dict = {
                "method": "cellsweep",
                "top_genes": parse_cell_filtering_table(files['cs_table'], "cellsweep"),
                "plot1": encode_image(files.get('cs_plot_contrib')),
                "plot2": encode_image(files.get('cs_plot_umap'))
            }
        elif files.get('cb_table'):
            cf_dict = {
                "method": "cellbender",
                "top_genes": parse_cell_filtering_table(files['cb_table'], "cellbender"),
                "plot1": encode_image(files.get('cb_plot1')),
                "plot2": encode_image(files.get('cb_plot2'))
            }

        if cf_dict: cell_filtering_data[s_id] = cf_dict

        # Merge Configs
        combined_config = raw_run_config.copy()
        if base_id in samplesheet_config:
            combined_config.update(samplesheet_config[base_id])

        samples_json_list.append({
            "sample_id": s_id,
            "config": combined_config,
            "source": "starsolo" if using_star else "alevin",
            "mapping": {
                "n_reads_sample": fmt(n_input_reads), "pct_r1_q30": pct_r1_q30, "pct_r2_q30": pct_r2_q30,
                "n_uniquely_mapped": fmt(n_unique), "pct_uniquely_mapped": pct_unique, "pct_multi_mapped": pct_multi,
                "pct_multi_too_many": pct_multi_too_many, "pct_unmapped_short": pct_short, "pct_unmapped_other": pct_other,
                "noise_pct": noise_pct, "intronic_pct": intronic_pct, "rrna_pct": rrna_pct,
                "mtdna_unique_pct": mtdna_unique, "mtdna_multi_all_pct": mtdna_multi_all, "mtdna_multi_primary_pct": mtdna_multi_primary,
                "n_cells": fmt(n_cells), "saturation": saturation, "reads_07_saturation": reads_07_sat_val,
                "mean_reads_per_cell": fmt(mean_reads), "median_umis": fmt(median_umis), "median_genes": fmt(median_genes),
                "total_genes": fmt(total_genes)
            },
            "cell_calling": {
                "expected_cells": fmt(samplesheet_config.get(base_id, {}).get("expected_cells", "N/A")),
                "num_cells": fmt(n_cells), "umi_threshold": umi_threshold
            },
            "taxonomy_sankey": extract_sankey_data(files.get('sankey'))
        })

        if files.get('per_cell'):
            try:
                with open(files['per_cell'], 'r') as f: per_cell_data[s_id] = json.load(f)
            except Exception: pass

        knee_data[s_id] = []
        if files.get('knee'):
            try:
                with open(files['knee'], 'r') as f:
                    knee_data[s_id] = [int(line.strip()) for line in f if line.strip().isdigit()]
            except Exception: pass

        saturation_images[s_id] = {}
        if files.get('sat_img'): saturation_images[s_id]["saturation"] = encode_image(files['sat_img'])
        if files.get('res_img'): saturation_images[s_id]["residuals"] = encode_image(files['res_img'])

    # Inject into HTML
    with open(args.template, 'r') as f: html_content = f.read()

    replacements = {
        "__RUN_METADATA_PLACEHOLDER__": json.dumps(report_metadata, indent=2),
        "__GLOBAL_DATA_PLACEHOLDER__": json.dumps({"overview": {"columns": global_cols, "rows": global_rows}}, indent=2),
        "__SAMPLES_DATA_PLACEHOLDER__": json.dumps({"samples": samples_json_list}, indent=2),
        "__PER_CELL_DATA_PLACEHOLDER__": json.dumps(per_cell_data, indent=2),
        "__SATURATION_IMAGES_PLACEHOLDER__": json.dumps(saturation_images, indent=2),
        "__KNEE_DATA_PLACEHOLDER__": json.dumps(knee_data, indent=2),
        "__CELLFILTERING_DATA_PLACEHOLDER__": json.dumps(cell_filtering_data, indent=2),
    }

    for placeholder, json_str in replacements.items():
        if placeholder in html_content:
            html_content = html_content.replace(placeholder, json_str)
        else:
            safe_ph = re.escape(placeholder)
            html_content = re.sub(r'\s*' + safe_ph + r'\s*', "\n" + json_str + "\n", html_content)

    with open(args.output, 'w') as f:
        f.write(html_content)

    print(f"Successfully generated {args.output} for {len(global_rows)} active analytical samples.", file=sys.stderr)

if __name__ == "__main__":
    main()
