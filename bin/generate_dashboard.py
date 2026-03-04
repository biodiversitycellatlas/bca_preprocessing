#!/usr/bin/env python3

import argparse
import json
import csv
import base64
import re
import os
import sys
import datetime

def parse_args():
    parser = argparse.ArgumentParser(description="Generate HTML Dashboard")
    parser.add_argument("--template", required=True, help="Path to HTML template")
    parser.add_argument("--output", required=True, help="Output HTML filename")

    # Core Metadata / Dashboard Headers
    parser.add_argument("--project", default="Biodiversity Cell Atlas", help="Project Name")
    parser.add_argument("--pipeline", default="bca-preprocessing", help="Pipeline Name")
    parser.add_argument("--version", default="0.2.0", help="Pipeline Version")
    parser.add_argument("--commit", default="unknown", help="Git Commit Hash")

    # Data Inputs
    parser.add_argument("--run_config", help="Path to run_config_date.txt")
    parser.add_argument("--samplesheet", required=True, help="Path to samplesheet.csv")

    # Data Sources
    parser.add_argument("--star_logs", nargs='*', help="List of STAR Log.final.out files")
    parser.add_argument("--star_summaries", nargs='*', help="List of STARsolo Summary.csv files")
    parser.add_argument("--star_full_logs", nargs='*', help="List of STAR Log.out files (for UMI threshold)")
    parser.add_argument("--mt_rrna_metrics", nargs='*', help="List of *_mt_rrna_metrics.txt files")
    parser.add_argument("--saturation_logs", nargs='*', help="List of saturation.log files")
    parser.add_argument("--cell_stats", nargs='*', help="List of STARsolo CellReads.stats files")

    # Visualization Files
    parser.add_argument("--sankey_files", nargs='*', help="List of *_kraken.sankey.html files")
    parser.add_argument("--per_cell_files", nargs='*', help="List of per-cell JSON files (e.g. *_metrics.json)")
    parser.add_argument("--saturation_imgs", nargs='*', help="List of saturation PNG images")
    parser.add_argument("--residuals_imgs", nargs='*', help="List of residuals PNG images")
    parser.add_argument("--knee_files", nargs='*', help="List of UMIperCellSorted.txt files")

    return parser.parse_args()

def find_file_for_sample(sample_id, file_list):
    """Finds a file in a list where the path matches the sample_id."""
    if not file_list:
        return None

    # 1. Exact filename match
    for f in file_list:
        filename = os.path.basename(f)
        if sample_id in filename:
            return f

    # 2. Directory match
    for f in file_list:
        if f"/{sample_id}/" in f or f"{os.sep}{sample_id}{os.sep}" in f:
            return f

    # 3. Loose match
    for f in file_list:
        if sample_id in f:
            return f

    return None

def to_pct(val, precision=2):
    """
    Converts a float (0.757819) to a percentage string (75.78%).
    Handles 'N/A' or None by returning them as-is.
    """
    if val is None or str(val).strip().upper() == "N/A" or val == "":
        return "N/A"
    try:
        float_val = float(val)
        return f"{float_val * 100:.{precision}f}%"
    except (ValueError, TypeError):
        return str(val)

def encode_image(image_path):
    if not image_path or not os.path.exists(image_path):
        return None
    try:
        with open(image_path, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
        return f"data:image/png;base64,{encoded_string}"
    except Exception as e:
        sys.stderr.write(f"Error encoding image {image_path}: {e}\n")
        return None

def extract_sankey_data(html_path):
    if not html_path or not os.path.exists(html_path):
        return None
    try:
        with open(html_path, 'r') as f:
            content = f.read()
            match = re.search(r'<script type="application/json" data-for="htmlwidget-[a-f0-9]+">(.*?)</script>', content, re.DOTALL)
            if match:
                data = json.loads(match.group(1))
                return data.get('x', None)
    except Exception as e:
        sys.stderr.write(f"Error extracting Sankey from {html_path}: {e}\n")
    return None


def parse_saturation_log(path: Optional[str]) -> str:
    """
    Extract read input needed to achieve n saturation from saturation.log.
    Returns a string like '22.5 M' or 'N/A'.
    """
    if not path or not os.path.exists(path):
        return "N/A"
    try:
        with open(path, 'r') as f:
            content = f.read()
        match = re.search(
            r"To achieve a saturation of [0-9.]+, ninput should be approximately:\s*(.*?)\s*reads",
            content
        )
        if match:
            return match.group(1).strip()
    except Exception as e:
        sys.stderr.write(f"Error parsing saturation log {path}: {e}\n")
    return "N/A"

def parse_run_config(path):
    config = {}
    if not path or not os.path.exists(path):
        return config
    try:
        with open(path, 'r') as f:
            for line in f:
                if '=' in line:
                    key, val = line.split('=', 1)
                    config[key.strip()] = val.strip()
    except Exception as e:
        sys.stderr.write(f"Error parsing run config: {e}\n")
    return config

def fmt(val):
    if val is None or val == "": return ""
    try:
        s_val = str(val).strip()
        if s_val.isdigit():
            return "{:,}".format(int(s_val))
    except:
        pass
    return str(val)

# --- Parsers ---

def parse_star_log(path):
    data = {}
    if not path or not os.path.exists(path):
        return data
    try:
        with open(path, 'r') as f:
            for line in f:
                if '|' in line:
                    parts = line.split('|')
                    key = parts[0].strip()
                    val = parts[1].strip()
                    data[key] = val
    except Exception as e:
        sys.stderr.write(f"Error parsing STAR log {path}: {e}\n")
    return data

def parse_starsolo_summary(path):
    data = {}
    if not path or not os.path.exists(path):
        return data
    try:
        with open(path, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) >= 2:
                    data[row[0].strip()] = row[1].strip()
    except Exception as e:
        sys.stderr.write(f"Error parsing STARsolo summary {path}: {e}\n")
    return data

def parse_mt_rrna_metrics(path):
    data = {}
    if not path or not os.path.exists(path):
        return data
    try:
        with open(path, 'r') as f:
            for line in f:
                if ',' in line:
                    parts = line.strip().split(',', 1)
                    data[parts[0].strip()] = parts[1].strip()
    except Exception as e:
        sys.stderr.write(f"Error parsing MT/rRNA metrics {path}: {e}\n")
    return data

def parse_starsolo_intronic(path):
    """
    Parses CellReads.stats to calculate intronic percentage of mapped reads.
    Formula: (nIntronic + nExonicIntronic) / (nMultimapped + nExonic + nIntronic + nExonicIntronic)
    """
    if not path or not os.path.exists(path):
        return "N/A"

    try:
        mapped_reads = 0
        intronic_sum = 0

        with open(path, 'r') as f:
            next(f) # Skip header
            for line in f:
                parts = line.strip().split()
                if len(parts) < 6:
                    continue

                # Column index: 2=Multi, 3=Exonic, 4=Intronic, 5=ExInt
                multi    = int(parts[2])
                exonic   = int(parts[3])
                intronic = int(parts[4])
                ex_int   = int(parts[5])

                intronic_sum += (intronic + ex_int)
                mapped_reads += (multi + exonic + intronic + ex_int)

        if mapped_reads > 0:
            return intronic_sum / mapped_reads
    except Exception as e:
        sys.stderr.write(f"Error parsing CellReads.stats {path}: {e}\n")

    return "N/A"

def extract_umi_cutoff(log_out_path, cfg_name="GeneFull_Ex50pAS"):
    if not log_out_path or not os.path.exists(log_out_path):
        return 0
    try:
        with open(log_out_path, 'r') as f:
            content = f.read()

        in_block = False
        for line in content.splitlines():
            if f"Starting Solo post-map for {cfg_name}" in line:
                in_block = True
                continue
            if in_block and "Starting Solo post-map for" in line:
                break
            if in_block and "cellFiltering" in line:
                match = re.search(r"nUMImin=(\d+)", line)
                if match:
                    return int(match.group(1))
    except Exception as e:
        sys.stderr.write(f"[WARNING] Failed to parse UMI cutoff from {log_out_path}: {e}\n")
    return 0

def main():
    args = parse_args()

    # 1. Load Raw Runtime Config (Files, Paths)
    raw_run_config = parse_run_config(args.run_config)

    # 2. Construct Clean Report Metadata
    report_metadata = {
        "project": args.project,
        "pipeline": args.pipeline,
        "version": args.version,
        "date": datetime.date.today().isoformat(),
        "commit": args.commit
    }

    # 3. Parse Samplesheet
    samplesheet_config = {}
    sample_ids = []

    if args.samplesheet and os.path.exists(args.samplesheet):
        with open(args.samplesheet, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get('sample'):
                    s_id = row['sample']
                    samplesheet_config[s_id] = row
                    sample_ids.append(s_id)
    else:
        sys.exit("Error: Samplesheet not found.")

    # 4. Process Samples
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

    print(f"Loaded {len(sample_ids)} samples from samplesheet.", file=sys.stderr)

    for s_id in sample_ids:
        # Find Core Files
        log_file = find_file_for_sample(s_id, args.star_logs)
        sum_file = find_file_for_sample(s_id, args.star_summaries)

        if not log_file or not sum_file:
            print(f"Skipping {s_id}: STAR log or Summary not found.", file=sys.stderr)
            continue

        print(f"Processing {s_id}...", file=sys.stderr)

        # Parse Basic Metrics
        star_stats = parse_star_log(log_file)
        solo_stats = parse_starsolo_summary(sum_file)

        mt_file = find_file_for_sample(s_id, args.mt_rrna_metrics)
        mt_stats = parse_mt_rrna_metrics(mt_file)

        # Parse Saturation Log
        sat_log_file = find_file_for_sample(s_id, args.saturation_logs)
        reads_07_sat_val = parse_saturation_log(sat_log_file)

        # Parse Intronic Metric from CellReads.stats
        cell_stats_file = find_file_for_sample(s_id, args.cell_stats)
        raw_intronic = parse_starsolo_intronic(cell_stats_file)
        intronic_pct = to_pct(raw_intronic)

        # UMI Threshold Detection
        umi_threshold = 0
        full_log_file = find_file_for_sample(s_id, args.star_full_logs)

        if not full_log_file and log_file:
            # Attempt to guess location
            guess_1 = log_file.replace("Log.final.out", "Log.out")
            guess_2 = log_file.replace("Log.final.out", f"{s_id}_Log.out")
            if os.path.exists(guess_1): full_log_file = guess_1
            elif os.path.exists(guess_2): full_log_file = guess_2

        if full_log_file and os.path.exists(full_log_file):
            umi_threshold = extract_umi_cutoff(full_log_file, cfg_name="GeneFull_Ex50pAS")

        # Helpers
        def get_val(source, key, default="N/A"):
            return source.get(key, default)

        # Metrics Extraction
        n_input_reads = get_val(star_stats, "Number of input reads")
        n_unique = get_val(star_stats, "Uniquely mapped reads number")
        pct_unique = get_val(star_stats, "Uniquely mapped reads %")
        pct_multi = get_val(star_stats, "Number of reads mapped to multiple loci")

        # STAR output format variance handling
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
        try:
            if frac_in_cells != "N/A":
                # Convert fraction to Noise %
                f_val = float(frac_in_cells)
                noise_pct = f"{(1.0 - f_val) * 100:.2f}%"
        except:
            pass

        # MT / rRNA Metrics
        rrna_pct = to_pct(get_val(mt_stats, "Percentage of rRNA reads (of uniquely mapped reads)"))
        mtdna_unique = to_pct(get_val(mt_stats, "Percentage of mtDNA reads (of mapped reads)"))
        mtdna_multi_all = to_pct(get_val(mt_stats, "Percentage of mtDNA in multimapped reads (all alignments)"))
        mtdna_multi_primary = to_pct(get_val(mt_stats, "Percentage of mtDNA in multimapped reads (primary alignment)"))

        # Global Row Construction
        g_row = [
            s_id,
            pct_unique,
            fmt(n_cells),
            saturation,
            reads_07_sat_val,
            noise_pct,
            fmt(median_umis),
            intronic_pct,
            rrna_pct,
            mtdna_unique,
            mtdna_multi_all
        ]
        global_rows.append(g_row)

        # Optional Files
        sankey_file = find_file_for_sample(s_id, args.sankey_files)
        per_cell_file = find_file_for_sample(s_id, args.per_cell_files)
        sat_img = find_file_for_sample(s_id, args.saturation_imgs)
        res_img = find_file_for_sample(s_id, args.residuals_imgs)
        knee_file = find_file_for_sample(s_id, args.knee_files)

        # JSON Object for Sample Detail View
        combined_config = raw_run_config.copy()
        if s_id in samplesheet_config:
            combined_config.update(samplesheet_config[s_id])

        sample_obj = {
            "sample_id": s_id,
            "config": combined_config,
            "mapping": {
                "n_reads_sample":        fmt(n_input_reads),
                "pct_r1_q30":            pct_r1_q30,
                "pct_r2_q30":            pct_r2_q30,
                "n_uniquely_mapped":     fmt(n_unique),
                "pct_uniquely_mapped":   pct_unique,
                "pct_multi_mapped":      pct_multi,
                "pct_multi_too_many":    pct_multi_too_many,
                "pct_unmapped_short":    pct_short,
                "pct_unmapped_other":    pct_other,
                "noise_pct":             noise_pct,
                "intronic_pct":          intronic_pct,
                "rrna_pct":              rrna_pct,
                "mtdna_unique_pct":      mtdna_unique,
                "mtdna_multi_all_pct":   mtdna_multi_all,
                "mtdna_multi_primary_pct": mtdna_multi_primary,
                "n_cells":               fmt(n_cells),
                "saturation":            saturation,
                "reads_07_saturation":   reads_07_sat_val,
                "mean_reads_per_cell":   fmt(mean_reads),
                "median_umis":           fmt(median_umis),
                "median_genes":          fmt(median_genes),
                "total_genes":           fmt(total_genes)
            },
            "cell_calling": {
                "expected_cells": fmt(samplesheet_config[s_id].get("expected_cells", "N/A")),
                "num_cells":      fmt(n_cells),
                "umi_threshold":  umi_threshold
            },
            "taxonomy_sankey": extract_sankey_data(sankey_file) if sankey_file else None
        }
        samples_json_list.append(sample_obj)

        # Store Large Data - Per Cell
        if per_cell_file:
            try:
                with open(per_cell_file, 'r') as f:
                    per_cell_data[s_id] = json.load(f)
            except Exception as e:
                print(f"  -> Error loading per-cell JSON for {s_id}: {e}", file=sys.stderr)
                per_cell_data[s_id] = {}
        else:
            per_cell_data[s_id] = {}

        # Knee Data
        if knee_file:
            try:
                with open(knee_file, 'r') as f:
                    knee_data[s_id] = [int(line.strip()) for line in f if line.strip().isdigit()]
            except:
                knee_data[s_id] = []

        # Images
        saturation_images[s_id] = {}
        if sat_img: saturation_images[s_id]["saturation"] = encode_image(sat_img)
        if res_img: saturation_images[s_id]["residuals"] = encode_image(res_img)

    # 5. Final Data Structure Assembly

    # Structure A: Global Overview
    global_data_struct = {
        "overview": {
            "columns": global_cols,
            "rows": global_rows
        }
    }

    # Structure B: Sample Details
    samples_data_struct = {"samples": samples_json_list}

    # 6. Read and Inject into Template
    with open(args.template, 'r') as f:
        html_content = f.read()

    replacements = {
        "__RUN_METADATA_PLACEHOLDER__": json.dumps(report_metadata, indent=2),
        "__GLOBAL_DATA_PLACEHOLDER__": json.dumps(global_data_struct, indent=2),
        "__SAMPLES_DATA_PLACEHOLDER__": json.dumps(samples_data_struct, indent=2),
        "__PER_CELL_DATA_PLACEHOLDER__": json.dumps(per_cell_data, indent=2),
        "__SATURATION_IMAGES_PLACEHOLDER__": json.dumps(saturation_images, indent=2),
        "__KNEE_DATA_PLACEHOLDER__": json.dumps(knee_data, indent=2),
    }

    for placeholder, json_str in replacements.items():
        if placeholder in html_content:
            html_content = html_content.replace(placeholder, json_str)
        else:
            # Fallback regex if spacing varies around the placeholder
            safe_ph = re.escape(placeholder)
            pattern = re.compile(r'\s*' + safe_ph + r'\s*')
            html_content = pattern.sub(lambda m: "\n" + json_str + "\n", html_content)

    with open(args.output, 'w') as f:
        f.write(html_content)

    print(f"Successfully generated {args.output} for {len(global_rows)} samples.", file=sys.stderr)

if __name__ == "__main__":
    main()
