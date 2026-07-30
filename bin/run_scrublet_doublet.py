#!/usr/bin/env python3
"""
Wraps Demuxafy's Scrublet.py CLI wrapper to add a fallback doublet
threshold: if Scrublet's automatic bimodal-threshold detection fails
(flags 0% or 100% of cells as doublets -- common on non-bimodal
datasets), Scrublet.py is re-run with a manual threshold (-t) set to the
90th percentile of the first pass's doublet scores.
"""

import argparse
import csv
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("RunScrublet")

FALLBACK_PERCENTILE = 90
BARCODE_COL = "Barcode"
CLASS_COL = "scrublet_DropletType"
SCORE_COL = "scrublet_Scores"


def run_scrublet(counts_matrix, outdir, threshold=None, filtered_barcodes=None):
    cmd = ["Scrublet.py", "-m", str(counts_matrix), "-o", str(outdir)]
    if threshold is not None:
        cmd += ["-t", str(threshold)]
    if filtered_barcodes is not None:
        cmd += ["-f", str(filtered_barcodes)]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def read_results(results_path):
    with open(results_path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def restrict_to_filtered_barcodes(results_path, filtered_barcodes_path):
    """
    Scrublet.py's own output keeps every original barcode (mostly 'unclassified'
    placeholders for the ones excluded by -f), even though only the provided
    candidate barcodes were actually scored. Prune the results file down to just
    those candidates so it's scoped the same way as scDblFinder's output (already
    restricted to just the filtered barcodes) and downstream consensus logic isn't
    swamped by tens of thousands of never-evaluated empty-droplet placeholder rows.
    """
    if filtered_barcodes_path is None:
        return

    with open(filtered_barcodes_path) as f:
        keep = {line.strip() for line in f if line.strip()}

    rows = read_results(results_path)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    kept_rows = [row for row in rows if row[BARCODE_COL] in keep]

    with open(results_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept_rows)

    logger.info(f"Restricted scrublet_results.tsv to {len(kept_rows)}/{len(rows)} provided candidate barcodes")


def doublet_fraction(rows):
    if not rows:
        return 0.0
    n_doublets = sum(1 for row in rows if row[CLASS_COL].strip().lower() == "doublet")
    return n_doublets / len(rows)


def parse_scores(rows):
    # Scrublet can leave the score blank for cells it couldn't score; skip those rather
    # than letting float('') crash the percentile fallback below.
    return np.array([float(row[SCORE_COL]) for row in rows if row[SCORE_COL].strip() != ""])


def main():
    parser = argparse.ArgumentParser(
        description="Run Demuxafy's Scrublet.py with a doublet-score-percentile fallback threshold."
    )
    parser.add_argument("--counts_matrix", type=str, required=True, help="10x-format matrix directory")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument("--sample_id", type=str, required=True, help="Sample ID (for logging)")
    parser.add_argument("--filtered_barcodes", type=str, default=None,
                         help="Optional list of real-cell barcodes to restrict doublet calling to "
                              "(passed through to Scrublet.py -f)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    run_scrublet(args.counts_matrix, outdir, filtered_barcodes=args.filtered_barcodes)

    results_path = outdir / "scrublet_results.tsv"
    restrict_to_filtered_barcodes(results_path, args.filtered_barcodes)
    rows = read_results(results_path)
    fraction = doublet_fraction(rows)

    if fraction in (0.0, 1.0):
        scores = parse_scores(rows)
        if scores.size == 0:
            logger.warning(
                f"[{args.sample_id}] Scrublet's automatic threshold flagged "
                f"{fraction * 100:.0f}% of cells as doublets, but no cells have a usable "
                f"doublet score to compute a fallback threshold from. Keeping the original result."
            )
        else:
            fallback_threshold = float(np.percentile(scores, FALLBACK_PERCENTILE))
            logger.warning(
                f"[{args.sample_id}] Scrublet's automatic threshold flagged "
                f"{fraction * 100:.0f}% of cells as doublets (likely a non-bimodal "
                f"score distribution). Re-running with a fallback threshold set to the "
                f"{FALLBACK_PERCENTILE}th percentile of doublet scores: {fallback_threshold:.4f}"
            )
            run_scrublet(args.counts_matrix, outdir, threshold=fallback_threshold, filtered_barcodes=args.filtered_barcodes)
            restrict_to_filtered_barcodes(results_path, args.filtered_barcodes)
            rows = read_results(results_path)
            fraction = doublet_fraction(rows)

    logger.info(f"[{args.sample_id}] Final doublet fraction: {fraction * 100:.2f}% ({len(rows)} cells evaluated)")


if __name__ == "__main__":
    main()
