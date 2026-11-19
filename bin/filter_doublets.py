#!/usr/bin/env python3
"""
Drops cells flagged as consensus doublets (the Scrublet and scDblFinder calls
combined by demuxafy_combine_results.R, under the --method given as
params.doublet_consensus_method; one tool's calls alone when the other could not
run) from a count-matrix triplet, writing a new triplet in its place.

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


# Display name and colour of each tool, by its column prefix in the combined table
TOOL_STYLE = {"scrublet": ("Scrublet", "#1f77b4"), "scDblFinder": ("scDblFinder", "#ff7f0e")}


def plot_filtering_summary(combined, doublets, method, n_evaluated, n_doublets, image_prefix):
    """
    Two-panel summary: which tool(s) called each doublet and which of those `method`
    removed, and the resulting before/after cell counts among the evaluated candidate cells.
    With one tool's calls alone (the other could not run), panel A shows that tool only.
    """
    present = doublet_calls.callers_present(combined)
    called = {tool: combined[f"{tool}_DropletType"].astype(str).str.lower() == "doublet" for tool in present}
    n_remaining = n_evaluated - n_doublets
    removed = combined["Barcode"].isin(doublets)

    # Constrained layout, not tight_layout: it accounts for panel A's equal aspect,
    # which tight_layout does not, and that clipped panel A's two-line title.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), layout="constrained")

    if len(present) == 2:
        scrublet_doublet, scdblfinder_doublet = called["scrublet"], called["scDblFinder"]
        both = int((scrublet_doublet & scdblfinder_doublet).sum())
        scrublet_only = int((scrublet_doublet & ~scdblfinder_doublet).sum())
        scdblfinder_only = int((~scrublet_doublet & scdblfinder_doublet).sum())

        # Name the removed set from the calls themselves rather than from the method's
        # rule, so methods beyond AnySinglet (intersection) and AnyDoublet (union) still
        # get a truthful title. Intersection first: when both sets coincide it reads best.
        if removed.equals(scrublet_doublet & scdblfinder_doublet):
            removed_label = f"overlap = removed by {method}"
        elif removed.equals(scrublet_doublet | scdblfinder_doublet):
            removed_label = f"union = removed by {method}"
        else:
            removed_label = f"{int(removed.sum())} removed by {method}"

        # Panel A: schematic (not area-proportional) overlap of per-method doublet calls
        ax1.add_patch(Circle((0.35, 0.5), 0.3, alpha=0.4, color=TOOL_STYLE["scrublet"][1]))
        ax1.add_patch(Circle((0.65, 0.5), 0.3, alpha=0.4, color=TOOL_STYLE["scDblFinder"][1]))
        ax1.text(0.22, 0.5, str(scrublet_only), ha="center", va="center", fontsize=13, fontweight="bold")
        ax1.text(0.78, 0.5, str(scdblfinder_only), ha="center", va="center", fontsize=13, fontweight="bold")
        ax1.text(0.5, 0.5, str(both), ha="center", va="center", fontsize=13, fontweight="bold")
        ax1.text(0.35, 0.88, TOOL_STYLE["scrublet"][0], ha="center", fontsize=11, color=TOOL_STYLE["scrublet"][1])
        ax1.text(0.65, 0.88, TOOL_STYLE["scDblFinder"][0], ha="center", fontsize=11, color=TOOL_STYLE["scDblFinder"][1])
        ax1.set_title(f"Doublets called per tool\n({removed_label})")
    else:
        # Panel A: the one tool that ran; every rule reduces to its own calls
        tool = present[0]
        name, colour = TOOL_STYLE[tool]
        missing = next(TOOL_STYLE[t][0] for t in TOOL_STYLE if t != tool)
        ax1.add_patch(Circle((0.5, 0.5), 0.3, alpha=0.4, color=colour))
        ax1.text(0.5, 0.5, str(int(called[tool].sum())), ha="center", va="center", fontsize=13, fontweight="bold")
        ax1.text(0.5, 0.88, name, ha="center", fontsize=11, color=colour)
        ax1.set_title(f"Doublets called by {name} alone\n({missing} could not run; removed by {method})")
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_aspect("equal")
    ax1.axis("off")

    # Panel B: before vs after cell counts among the evaluated candidate cells
    counts = [n_evaluated, n_remaining]
    bars = ax2.bar(["Before filtering", "After filtering"], counts, color=["#aec7e8", "#2ca02c"])
    for bar, value in zip(bars, counts):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(value), ha="center", va="bottom")
    ax2.set_ylabel("Number of cells")
    pct = (n_doublets / n_evaluated * 100) if n_evaluated > 0 else 0.0
    ax2.set_title(f"{n_doublets} {method} doublets removed ({pct:.1f}%)")

    plt.savefig(f"{image_prefix}doublet_filtering_summary.png", dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Filter consensus doublets out of a count-matrix triplet.")
    parser.add_argument("--mtx", type=str, required=True, help="Matrix to filter")
    parser.add_argument("--barcodes", type=str, required=True, help="Barcode axis of --mtx")
    parser.add_argument("--features", type=str, required=True, help="Feature axis of --mtx")
    parser.add_argument("--combined_results", type=str, required=True,
                         help="demuxafy_combine_results.R '_w_combined_assignments.tsv' output")
    parser.add_argument("--method", type=str, required=True,
                         help="Consensus method name used in demuxafy_combine_results.R (e.g. AnySinglet)")
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

    plot_filtering_summary(combined, doublets, args.method, n_evaluated, n_doublets, args.image_prefix)

    logger.info("Doublet filtering complete.")


if __name__ == "__main__":
    main()
