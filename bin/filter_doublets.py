#!/usr/bin/env python3
"""
Drops cells flagged as consensus doublets (Scrublet + scDblFinder agreement,
combined via Demuxafy's Combine_Results.R) from a count-matrix triplet, writing a
new triplet in its place.

Works on the matrix rather than on an .h5ad because it sits in the middle of the
chain: the triplet it emits is what CellSweep denoises and what MTX_TO_H5AD
finally packs, so the published AnnData is the doublet-free matrix. Runs once on
the raw matrix and once on the cell-called one, and only when
params.perform_doublet_filtering is set -- otherwise the doublets are kept and
merely annotated (by CellSweep on the UMAP, and by MTX_TO_H5AD in .obs).
"""

import argparse
import logging
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import doublet_calls
import mtx_io

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
    parser = argparse.ArgumentParser(description="Filter consensus doublets out of a count-matrix triplet.")
    parser.add_argument("--mtx", type=str, required=True, help="Matrix to filter")
    parser.add_argument("--barcodes", type=str, required=True, help="Barcode axis of --mtx")
    parser.add_argument("--features", type=str, required=True, help="Feature axis of --mtx")
    parser.add_argument("--combined_results", type=str, required=True,
                         help="Combine_Results.R '_w_combined_assignments.tsv' output")
    parser.add_argument("--method", type=str, required=True,
                         help="Consensus method name used in Combine_Results.R (e.g. AnySinglet)")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory for the filtered triplet")
    parser.add_argument("--summary_txt", type=str, required=True, help="Output path for summary text file")
    parser.add_argument("--image_prefix", type=str, default="", help="Prefix for output PNGs")
    args = parser.parse_args()

    # CSC (Compressed Sparse Column) for efficient cell slicing: the matrix arrives
    # genes x cells whichever mapper it came from, so cells are its columns.
    matrix, barcodes, features = mtx_io.read_triplet(args.mtx, args.barcodes, args.features)
    matrix = matrix.tocsc()

    logger.info(f"Loading combined doublet results from {args.combined_results}")
    combined = doublet_calls.read_combined_results(args.combined_results)
    doublets = doublet_calls.doublet_barcodes(combined, args.method, source=args.combined_results)

    # On the raw matrix `barcodes` holds every droplet (incl. empty ones) while only the
    # cell-called ones in `combined` were evaluated; on the cell-called matrix the two sets
    # coincide. Report both counts either way rather than assuming which matrix this is.
    total_cells = len(barcodes)
    n_evaluated = len(combined)
    keep_idx = [i for i, barcode in enumerate(barcodes) if barcode not in doublets]
    n_doublets = total_cells - len(keep_idx)
    doublet_pct = (n_doublets / n_evaluated * 100) if n_evaluated > 0 else 0.0

    logger.info(
        f"Removing {n_doublets} / {n_evaluated} evaluated candidate cells as consensus doublets ({doublet_pct:.2f}%)"
    )

    if not keep_idx:
        raise SystemExit(
            f"Error: every one of the {total_cells} barcodes in {args.barcodes} was called a "
            f"consensus doublet; there is no matrix left to write."
        )

    # Barcodes keep their source order, so the filtered triplet stays comparable to the
    # matrix it came from.
    mtx_io.write_triplet(
        args.outdir,
        matrix[:, keep_idx].tocsr(),
        [barcodes[i] for i in keep_idx],
        features,
    )

    with open(args.summary_txt, "w") as f:
        f.write(f"Total_Cells_In_Input_Matrix\t{total_cells}\n")
        f.write(f"Cells_Evaluated_For_Doublets\t{n_evaluated}\n")
        f.write(f"Consensus_Doublets_Removed\t{n_doublets}\n")
        f.write(f"Doublet_Percentage\t{doublet_pct:.2f}%\n")
        f.write(f"Cells_Remaining\t{len(keep_idx)}\n")

    plot_filtering_summary(combined, n_evaluated, n_doublets, args.image_prefix)

    logger.info("Doublet filtering complete.")


if __name__ == "__main__":
    main()
