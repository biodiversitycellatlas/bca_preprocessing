#!/usr/bin/env python3
"""
Drops cells flagged as consensus doublets (Scrublet + scDblFinder agreement,
combined via Demuxafy's Combine_Results.R) from a raw .h5ad file before it
is passed downstream to CellSweep.
"""

import argparse
import logging
import sys

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("FilterDoublets")


def plot_filtering_summary(combined, n_evaluated, n_doublets, image_prefix):
    """
    Two-panel summary: which method(s) called each doublet (overlap = consensus),
    and the resulting before/after cell counts among the evaluated candidate cells.
    """
    scrublet_doublet = combined["scrublet_DropletType"].str.lower() == "doublet"
    scdblfinder_doublet = combined["scDblFinder_DropletType"].str.lower() == "doublet"

    both = int((scrublet_doublet & scdblfinder_doublet).sum())
    scrublet_only = int((scrublet_doublet & ~scdblfinder_doublet).sum())
    scdblfinder_only = int((~scrublet_doublet & scdblfinder_doublet).sum())
    n_remaining = n_evaluated - n_doublets

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # Panel A: schematic (not area-proportional) overlap of per-method doublet calls
    ax1.add_patch(Circle((0.35, 0.5), 0.3, alpha=0.4, color="#1f77b4"))
    ax1.add_patch(Circle((0.65, 0.5), 0.3, alpha=0.4, color="#ff7f0e"))
    ax1.text(0.22, 0.5, str(scrublet_only), ha="center", va="center", fontsize=13, fontweight="bold")
    ax1.text(0.78, 0.5, str(scdblfinder_only), ha="center", va="center", fontsize=13, fontweight="bold")
    ax1.text(0.5, 0.5, str(both), ha="center", va="center", fontsize=13, fontweight="bold")
    ax1.text(0.35, 0.88, "Scrublet", ha="center", fontsize=11, color="#1f77b4")
    ax1.text(0.65, 0.88, "scDblFinder", ha="center", fontsize=11, color="#ff7f0e")
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_aspect("equal")
    ax1.axis("off")
    ax1.set_title("Doublets called per method\n(overlap = consensus, removed)")

    # Panel B: before vs after cell counts among the evaluated candidate cells
    counts = [n_evaluated, n_remaining]
    bars = ax2.bar(["Before filtering", "After filtering"], counts, color=["#aec7e8", "#2ca02c"])
    for bar, value in zip(bars, counts):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(value), ha="center", va="bottom")
    ax2.set_ylabel("Number of cells")
    pct = (n_doublets / n_evaluated * 100) if n_evaluated > 0 else 0.0
    ax2.set_title(f"{n_doublets} consensus doublets removed ({pct:.1f}%)")

    plt.tight_layout()
    plt.savefig(f"{image_prefix}doublet_filtering_summary.png", dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Filter consensus doublets out of an h5ad file.")
    parser.add_argument("--input_h5ad", type=str, required=True, help="Path to raw h5ad file")
    parser.add_argument("--combined_results", type=str, required=True,
                         help="Combine_Results.R '_w_combined_assignments.tsv' output")
    parser.add_argument("--method", type=str, required=True,
                         help="Consensus method name used in Combine_Results.R (e.g. AnySinglet)")
    parser.add_argument("--output_h5ad", type=str, required=True, help="Output path for filtered h5ad")
    parser.add_argument("--summary_txt", type=str, required=True, help="Output path for summary text file")
    parser.add_argument("--image_prefix", type=str, default="", help="Prefix for output PNGs")
    args = parser.parse_args()

    logger.info(f"Loading {args.input_h5ad}")
    adata = ad.read_h5ad(args.input_h5ad)

    logger.info(f"Loading combined doublet results from {args.combined_results}")
    combined = pd.read_csv(args.combined_results, sep="\t")

    classification_col = f"{args.method}_DropletType"
    if classification_col not in combined.columns:
        raise ValueError(
            f"Expected column '{classification_col}' not found in {args.combined_results}. "
            f"Available columns: {list(combined.columns)}"
        )

    doublet_barcodes = set(
        combined.loc[combined[classification_col].str.lower() == "doublet", "Barcode"]
    )

    # adata still holds the full raw/full matrix (incl. empty droplets); only the candidate
    # cells in `combined` (e.g. STARsolo's filtered barcode list) were actually evaluated.
    total_cells = adata.n_obs
    n_evaluated = len(combined)
    is_doublet = adata.obs_names.isin(doublet_barcodes)
    n_doublets = int(is_doublet.sum())
    doublet_pct = (n_doublets / n_evaluated * 100) if n_evaluated > 0 else 0.0

    logger.info(
        f"Removing {n_doublets} / {n_evaluated} evaluated candidate cells as consensus doublets ({doublet_pct:.2f}%)"
    )

    adata_filtered = adata[~is_doublet].copy()
    adata_filtered.write_h5ad(args.output_h5ad)

    with open(args.summary_txt, "w") as f:
        f.write(f"Total_Cells_In_Raw_Matrix\t{total_cells}\n")
        f.write(f"Cells_Evaluated_For_Doublets\t{n_evaluated}\n")
        f.write(f"Consensus_Doublets_Removed\t{n_doublets}\n")
        f.write(f"Doublet_Percentage\t{doublet_pct:.2f}%\n")
        f.write(f"Cells_Remaining\t{adata_filtered.n_obs}\n")

    plot_filtering_summary(combined, n_evaluated, n_doublets, args.image_prefix)

    logger.info("Doublet filtering complete.")


if __name__ == "__main__":
    main()
