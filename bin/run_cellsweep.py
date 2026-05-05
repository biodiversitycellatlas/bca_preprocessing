#!/usr/bin/env python3
"""
Handles empty droplet detection, automated cell typing via Leiden clustering,
and ambient RNA denoising using the Cellsweep tool.
"""

import argparse
import logging
import os
import sys
from pathlib import Path

import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from cellsweep import denoise_count_matrix
from cellsweep.utils import visualization_utils as cs_utils


# Configure Logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("CellSweep")


def detect_empty_droplets(adata, expected_cells, image_prefix):
    """
    Identify empty droplets using the knee-plot method.
    """
    logger.info("Running empty droplet detection...")

    knee_path = f"{image_prefix}knee_plot.png"
    umi_cutoff = cs_utils.knee_plot(
        adata,
        expected_cells=expected_cells,
        out_path=knee_path,
        show=False
    )

    adata = cs_utils.infer_empty_droplets(
        adata,
        method="threshold",
        umi_cutoff=umi_cutoff
    )
    return adata


def assign_clusters(adata, min_genes, max_mt_percent, remove_doublets, image_prefix):
    """
    Perform clustering to provide CellSweep with biological groups.
    When remove_doublets is False, doublets are identified but retained;
    their predictions are propagated back to adata.obs for downstream visualization.
    """
    logger.info("Assigning preliminary cell types via Leiden clustering...")

    adata_celltype = adata[~adata.obs["is_empty"]].copy()

    if min_genes is not None:
        sc.pp.filter_cells(adata_celltype, min_genes=min_genes)

    sc.pp.filter_genes(adata_celltype, min_cells=1)

    if max_mt_percent is not None:
        adata_celltype.var["mt"] = adata_celltype.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata_celltype, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        adata_celltype = adata_celltype[adata_celltype.obs.pct_counts_mt < max_mt_percent, :]

    # Always run Scrublet — either to remove or to annotate doublets for visualization
    logger.info("Running Scrublet doublet detection...")
    sc.pp.scrublet(adata_celltype)

    total_cells = len(adata_celltype)
    num_doublets = int(adata_celltype.obs["predicted_doublet"].sum())
    doublet_pct = (num_doublets / total_cells * 100) if total_cells > 0 else 0.0

    action = "Removed" if remove_doublets else "Identified (kept for visualization)"
    logger.info(f"Doublets {action}: {num_doublets} ({doublet_pct:.2f}%).")

    doublet_out_path = f"{image_prefix}doublet_summary.txt"
    with open(doublet_out_path, "w") as f:
        f.write(f"Total_Cells_Evaluated\t{total_cells}\n")
        f.write(f"Doublets_{'Removed' if remove_doublets else 'Identified'}\t{num_doublets}\n")
        f.write(f"Doublet_Percentage\t{doublet_pct:.2f}%\n")

    if remove_doublets:
        adata_celltype = adata_celltype[~adata_celltype.obs["predicted_doublet"]]
    else:
        # Propagate doublet predictions to the full adata for downstream UMAP highlighting
        adata.obs["predicted_doublet"] = (
            adata_celltype.obs["predicted_doublet"]
            .reindex(adata.obs.index)
            .fillna(False)
            .astype(bool)
        )
        adata.obs["doublet_score"] = (
            adata_celltype.obs["doublet_score"]
            .reindex(adata.obs.index)
            .fillna(0.0)
        )

    # Standard Scanpy preprocessing for clustering
    adata_celltype.layers["counts"] = adata_celltype.X.copy()
    sc.pp.normalize_total(adata_celltype, target_sum=1e4)
    sc.pp.log1p(adata_celltype)
    sc.pp.highly_variable_genes(adata_celltype, n_top_genes=2000)
    sc.tl.pca(adata_celltype, svd_solver="arpack", random_state=42)
    sc.pp.neighbors(adata_celltype, n_neighbors=15, n_pcs=50)

    sc.tl.leiden(adata_celltype, flavor="igraph", n_iterations=2, resolution=1.0, random_state=42)

    adata.obs["celltype"] = (
        adata_celltype.obs["leiden"]
        .cat.add_categories(["empty"])
        .reindex(adata.obs.index)
        .fillna("empty")
    )

    adata_cells_filtered = adata[
        (adata.obs["is_empty"] | adata.obs_names.isin(adata_celltype.obs_names))
    ].copy()

    return adata


def generate_visualizations(adata_cs, image_prefix):
    """
    Generate QC plots for the denoising results.
    """
    logger.info("Generating diagnostic visualizations...")

    cs_utils.plot_histogram(
        adata_cs,
        col="alpha_hat",
        out_path=f"{image_prefix}alpha_hat_dist.png",
        kind="cdf"
    )

    cs_utils.plot_histogram(
        adata_cs, adata_df="var", col="ambient_hat",
        out_path=f"{image_prefix}ambient_hat_histogram.png", kind="pdf", ylog=True
    )
    cs_utils.print_top_ambient_genes(adata_cs, top_n=20, out_path=f"{image_prefix}top_ambient_genes.csv")


def generate_noise_boxplot(adata_cs, image_prefix, label):
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=adata_cs.obs, x='celltype', y='alpha_hat', palette='viridis')
    plt.title('Ambient Fraction by Cluster')
    plt.savefig(f"{image_prefix}contamination_{label}.png")
    plt.close('all')


def compare_umaps(adata, image_prefix, highlight_doublets=False):
    """
    Manifold comparison (Before vs. After CellSweep denoising).
    When highlight_doublets=True, renders a 2x2 grid showing both cell-type
    labels and doublet predictions for the raw and denoised embeddings.
    """
    logger.info("Computing UMAPs for comparison...")
    cell_mask = (adata.obs['celltype'] != 'empty')
    ad_plot = adata[cell_mask].copy()

    # Raw Embedding (from raw layer)
    ad_raw = ad_plot.copy()
    ad_raw.X = ad_raw.layers['raw']
    sc.pp.normalize_total(ad_raw)
    sc.pp.log1p(ad_raw)
    sc.pp.pca(ad_raw)
    sc.pp.neighbors(ad_raw)
    sc.tl.umap(ad_raw)

    # Denoised Embedding (current X)
    sc.pp.pca(ad_plot)
    sc.pp.neighbors(ad_plot)
    sc.tl.umap(ad_plot)

    has_doublet_col = (
        highlight_doublets
        and "predicted_doublet" in ad_plot.obs.columns
    )

    if has_doublet_col:
        # Add a human-readable label column for clean legend text
        for ad_obj in (ad_raw, ad_plot):
            ad_obj.obs["doublet_label"] = (
                ad_obj.obs["predicted_doublet"]
                .map({True: "Doublet", False: "Singlet"})
                .astype("category")
            )

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))

        sc.pl.umap(
            ad_raw, color="celltype",
            title="Raw Counts – Cell Types",
            show=False, ax=axes[0, 0]
        )
        sc.pl.umap(
            ad_raw, color="doublet_label",
            title="Raw Counts – Predicted Doublets",
            palette={"Doublet": "#d62728", "Singlet": "#aec7e8"},
            show=False, ax=axes[0, 1]
        )
        sc.pl.umap(
            ad_plot, color="celltype",
            title="Denoised Counts – Cell Types",
            show=False, ax=axes[1, 0]
        )
        sc.pl.umap(
            ad_plot, color="doublet_label",
            title="Denoised Counts – Predicted Doublets",
            palette={"Doublet": "#d62728", "Singlet": "#aec7e8"},
            show=False, ax=axes[1, 1]
        )
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        sc.pl.umap(ad_raw, color='celltype', title='Raw Counts', show=False, ax=ax1)
        sc.pl.umap(ad_plot, color='celltype', title='Filtered Counts', show=False, ax=ax2)

    plt.tight_layout()
    plt.savefig(f"{image_prefix}umap_comparison.png", dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Run CellSweep denoising on scRNA-seq data.")
    parser.add_argument("--input_h5ad", type=str, required=True, help="Path to raw h5ad file")
    parser.add_argument("--expected_cells", type=int, required=True, help="Expected number of cells")
    parser.add_argument("--cs_filtered_h5ad", type=str, required=True, help="Output path for filtered cells")
    parser.add_argument("--cs_full_h5ad", type=str, required=True, help="Output path for full denoised matrix")
    parser.add_argument("--image_prefix", type=str, default="", help="Prefix for output PNGs")
    parser.add_argument("--alpha_threshold", type=float, default=0.3, help="Max ambient fraction allowed")
    parser.add_argument("--min-genes", type=int, default=10, help="Minimum number of genes")
    parser.add_argument("--max-mt-percent", type=int, default=25, help="Maximum percentage of MT content, scale 0 to 100.")
    parser.add_argument('--remove-doublets', default=True, action=argparse.BooleanOptionalAction,
                        help="Remove doublets before clustering (recommended). If disabled, doublets are "
                             "identified and highlighted on the UMAP instead.")
    parser.add_argument("--threads", type=int, default=1, help="Threads for CellSweep")
    args = parser.parse_args()

    # 1. Load Data
    logger.info(f"Loading data from {args.input_h5ad}")
    adata = sc.read_h5ad(args.input_h5ad)
    adata.var_names_make_unique()

    # 2. Identify Empty Droplets
    adata = detect_empty_droplets(adata, args.expected_cells, args.image_prefix)

    # 3. Cluster Cells (doublets always identified; optionally removed)
    adata = assign_clusters(adata, args.min_genes, args.max_mt_percent, args.remove_doublets, args.image_prefix)

    # 4. Run CellSweep Denoising
    logger.info("Starting CellSweep denoising (this may take a few minutes)...")
    log_file = f"{args.image_prefix}cellsweep_run.log"

    adata_cs = denoise_count_matrix(
        adata,
        adata_out=args.cs_full_h5ad,
        log_file=log_file,
        threads=args.threads
    )
    generate_noise_boxplot(adata_cs, args.image_prefix, 'before')

    # 5. Filter and Save
    generate_visualizations(adata_cs, args.image_prefix)

    # Pass highlight_doublets=True when doublets were kept instead of removed
    compare_umaps(adata_cs, args.image_prefix, highlight_doublets=not args.remove_doublets)

    logger.info(f"Filtering cells (is_empty=False, alpha_hat <= {args.alpha_threshold})")
    keep_mask = (~adata_cs.obs["is_empty"]) & (adata_cs.obs["alpha_hat"] <= args.alpha_threshold)
    adata_filtered = adata_cs[keep_mask].copy()
    generate_noise_boxplot(adata_filtered, args.image_prefix, 'after')

    adata_filtered.write_h5ad(args.cs_filtered_h5ad)
    logger.info("CellSweep complete. Output files written.")


if __name__ == "__main__":
    main()
