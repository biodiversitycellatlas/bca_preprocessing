#!/usr/bin/env python3
"""
Wraps demuxafy_scrublet.py (Demuxafy's Scrublet.py, ported) to add a fallback
doublet threshold: if Scrublet's automatic threshold lands at an extreme (flags
0% or 100% of cells as doublets), demuxafy_scrublet.py is re-run with a manual
threshold (-t) set to the 90th percentile of the first pass's doublet
scores.

That fallback assumes there are doublets to find and only the threshold
went astray, so it is applied only when Scrublet's own thresholding
succeeded, i.e. when it did find a bimodal split in the score
distribution.  When it reports no such split, forcing the 90th percentile
would label a tenth of the barcodes as doublets purely by construction,
so the 0% is reported as it stands instead.
"""

import argparse
import csv
import logging
import shutil
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
BOOTSTRAP_THRESHOLD = 1.0
BARCODE_COL = "Barcode"
CLASS_COL = "scrublet_DropletType"
SCORE_COL = "scrublet_Scores"

# demuxafy_scrublet.py writes "scrublet_results.tsv" at its automatic threshold and
# "scrublet_doublets_singlets.tsv" at a manual one (-t), as Demuxafy's Scrublet.py did.
# A manual re-run comes after the automatic one, so its file is looked for first; the
# result is handed on as "scrublet_results.tsv", which is what the Nextflow module collects.
RESULTS_CANDIDATES = ("scrublet_doublets_singlets.tsv", "scrublet_results.tsv")
CANONICAL_RESULTS = "scrublet_results.tsv"

# demuxafy_scrublet.py's exit status when its automatic threshold finds no bimodal split
NO_THRESHOLD_EXIT = 3


def run_scrublet(counts_matrix, outdir, doublet_rate_per_1000, threshold=None, filtered_barcodes=None):
    cmd = ["demuxafy_scrublet.py", "-m", str(counts_matrix), "-o", str(outdir),
           "--doublet_rate_per_1000", str(doublet_rate_per_1000)]
    if threshold is not None:
        cmd += ["-t", str(threshold)]
    if filtered_barcodes is not None:
        cmd += ["-f", str(filtered_barcodes)]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def resolve_results(outdir):
    """
    Locate the per-barcode results demuxafy_scrublet.py just wrote.

    Re-resolved after every call rather than cached, since a re-run at a manual
    threshold writes under the other name.
    """
    for name in RESULTS_CANDIDATES:
        candidate = outdir / name
        if candidate.exists():
            return candidate

    present = sorted(p.name for p in outdir.iterdir()) if outdir.is_dir() else []
    raise SystemExit(
        f"demuxafy_scrublet.py wrote none of {', '.join(RESULTS_CANDIDATES)} into {outdir}. "
        f"Found: {', '.join(present) if present else '(nothing)'}"
    )


def read_results(results_path):
    with open(results_path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def restrict_to_filtered_barcodes(results_path, filtered_barcodes_path):
    """
    Keep only the candidate barcodes in the results file, so it is scoped the same
    way as scDblFinder's output. demuxafy_scrublet.py already writes just those
    under -f; this guards the consensus against any never-evaluated rows.
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
        description="Run demuxafy_scrublet.py with a doublet-score-percentile fallback threshold."
    )
    parser.add_argument("--counts_matrix", type=str, required=True, help="10x-format matrix directory")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument("--sample_id", type=str, required=True, help="Sample ID (for logging)")
    parser.add_argument("--doublet_rate_per_1000", type=float, default=0.008,
                         help="Expected doublet rate per 1,000 cells (passed through to demuxafy_scrublet.py)")
    parser.add_argument("--filtered_barcodes", type=str, default=None,
                         help="Optional list of real-cell barcodes to restrict doublet calling to "
                              "(passed through to demuxafy_scrublet.py -f)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rate = args.doublet_rate_per_1000

    # Scrublet's own automatic thresholding doubles as the bimodality test: it is what
    # fails when the score distribution has no separable doublet mode.
    bimodal = True
    try:
        run_scrublet(args.counts_matrix, outdir, rate, filtered_barcodes=args.filtered_barcodes)
    except subprocess.CalledProcessError as err:
        # Any other failure is a real error: the task fails, and the consensus goes on
        # with scDblFinder's calls alone
        if err.returncode != NO_THRESHOLD_EXIT:
            raise
        # No bimodal split: Scrublet set no threshold and called nothing. Re-run with an
        # explicit threshold so the scores get written; that leaves 0% doublets, which
        # is then reported as it stands.
        bimodal = False
        logger.warning(
            f"[{args.sample_id}] Scrublet's automatic threshold detection found no bimodal "
            f"split. Re-running with a bootstrap threshold of {BOOTSTRAP_THRESHOLD} to "
            f"recover the doublet scores."
        )
        run_scrublet(args.counts_matrix, outdir, rate, threshold=BOOTSTRAP_THRESHOLD,
                     filtered_barcodes=args.filtered_barcodes)

    results_path = resolve_results(outdir)
    restrict_to_filtered_barcodes(results_path, args.filtered_barcodes)
    rows = read_results(results_path)
    fraction = doublet_fraction(rows)

    if fraction in (0.0, 1.0):
        scores = parse_scores(rows)
        if not bimodal:
            # Re-thresholding at a percentile of a unimodal score distribution does not
            # recover doublets that thresholding missed, it manufactures them: the 90th
            # percentile calls 10% of barcodes doublets whatever the scores look like.
            logger.warning(
                f"[{args.sample_id}] Scrublet found no bimodal split in the doublet score "
                f"distribution, so there is no threshold that separates doublets from "
                f"singlets here. Reporting the {fraction * 100:.0f}% it called rather than "
                f"forcing the {FALLBACK_PERCENTILE}th percentile, which would label "
                f"{100 - FALLBACK_PERCENTILE}% of barcodes as doublets by construction."
            )
        elif scores.size == 0:
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
            run_scrublet(args.counts_matrix, outdir, rate, threshold=fallback_threshold,
                         filtered_barcodes=args.filtered_barcodes)
            results_path = resolve_results(outdir)
            restrict_to_filtered_barcodes(results_path, args.filtered_barcodes)
            rows = read_results(results_path)
            fraction = doublet_fraction(rows)

    # Hand the results on under one name, whichever the image happened to write
    canonical_path = outdir / CANONICAL_RESULTS
    if results_path != canonical_path:
        shutil.copyfile(results_path, canonical_path)
        logger.info(f"Copied {results_path.name} to {CANONICAL_RESULTS}")

    logger.info(f"[{args.sample_id}] Final doublet fraction: {fraction * 100:.2f}% ({len(rows)} cells evaluated)")


if __name__ == "__main__":
    main()
