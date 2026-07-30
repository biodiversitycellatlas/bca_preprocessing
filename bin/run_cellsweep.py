#!/usr/bin/env python3
"""
Handles empty droplet detection, automated cell typing via Leiden clustering,
and ambient RNA denoising using the Cellsweep tool. Doublet detection happens
upstream on the cell-called matrix (Scrublet + scDblFinder consensus); the calls
are passed in here purely as an annotation to project onto the UMAP. They are
only passed when those cells are still present: removing them instead is opt-in
and handled upstream (see bin/filter_doublets.py), and in that case this script
runs without --doublet_results.
"""

import argparse
import logging
import sys

import numpy as np
import pandas as pd
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


DOUBLET_COL = "doublet_status"


def read_doublet_barcodes(combined_results, method):
    """
    Barcodes called doublets by the upstream consensus (Demuxafy's Combine_Results.R
    output for `method`, e.g. AnySinglet = the intersection of the per-tool calls).
    """
    combined = pd.read_csv(combined_results, sep="\t")

    classification_col = f"{method}_DropletType"
    if classification_col not in combined.columns:
        raise ValueError(
            f"Expected column '{classification_col}' not found in {combined_results}. "
            f"Available columns: {list(combined.columns)}"
        )

    return set(combined.loc[combined[classification_col].str.lower() == "doublet", "Barcode"])


def annotate_doublets(adata, doublet_barcodes):
    """
    Flag the consensus doublets in .obs. Doublets were called on the cell-called matrix,
    so only a subset of this (raw) matrix's barcodes was ever evaluated; everything else
    is left as 'singlet'. Applied again after denoising, since CellSweep hands back its
    own AnnData object.
    """
    is_doublet = adata.obs_names.isin(doublet_barcodes)
    adata.obs[DOUBLET_COL] = pd.Categorical(
        np.where(is_doublet, "doublet", "singlet"),
        categories=["singlet", "doublet"]
    )
    logger.info(f"Annotated {int(is_doublet.sum())} / {adata.n_obs} barcodes as consensus doublets")
    return adata


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


def assign_clusters(adata, min_genes, max_mt_percent):
    """
    Perform clustering to provide CellSweep with biological groups.
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


def compare_umaps(adata, image_prefix):
    """
    Manifold comparison (Raw vs. CellSweep-denoised counts), plus a projection of the
    consensus doublets onto the denoised embedding when they were annotated.
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

    has_doublets = DOUBLET_COL in ad_plot.obs
    n_panels = 3 if has_doublets else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 6))
    sc.pl.umap(ad_raw, color='celltype', title='Raw Counts', show=False, ax=axes[0])
    sc.pl.umap(ad_plot, color='celltype', title='Denoised Counts', show=False, ax=axes[1])

    if has_doublets:
        n_doublets = int((ad_plot.obs[DOUBLET_COL] == 'doublet').sum())
        # `groups` greys out the singlets and draws the doublets on top of them, but needs
        # a non-empty group -- fall back to a plain two-colour scatter when none are left
        # (e.g. when they were already dropped by params.perform_doublet_filtering).
        sc.pl.umap(
            ad_plot,
            color=DOUBLET_COL,
            groups=['doublet'] if n_doublets > 0 else None,
            palette={'singlet': 'lightgrey', 'doublet': '#d62728'},
            title=f'Consensus doublets (n={n_doublets})',
            show=False,
            ax=axes[2]
        )

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
    parser.add_argument("--threads", type=int, default=1, help="Threads for CellSweep")
    parser.add_argument("--doublet_results", type=str, default=None,
                        help="Optional Combine_Results.R '_w_combined_assignments.tsv' output; "
                             "its consensus doublets are annotated and projected onto the UMAP")
    parser.add_argument("--doublet_method", type=str, default=None,
                        help="Consensus method name used in Combine_Results.R (e.g. AnySinglet). "
                             "Required with --doublet_results")
    args = parser.parse_args()

    if args.doublet_results and not args.doublet_method:
        parser.error("--doublet_method is required when --doublet_results is given")

    # 1. Load Data
    logger.info(f"Loading data from {args.input_h5ad}")
    adata = sc.read_h5ad(args.input_h5ad)
    adata.var_names_make_unique()

    # 1b. Annotate (do not remove) the consensus doublets called upstream
    doublet_barcodes = None
    if args.doublet_results:
        logger.info(f"Loading consensus doublet calls from {args.doublet_results}")
        doublet_barcodes = read_doublet_barcodes(args.doublet_results, args.doublet_method)
        adata = annotate_doublets(adata, doublet_barcodes)

    # 2. Identify Empty Droplets
    adata = detect_empty_droplets(adata, args.expected_cells, args.image_prefix)

    # 3. Cluster Cells (biological grouping for CellSweep)
    adata = assign_clusters(adata, args.min_genes, args.max_mt_percent)

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

    if doublet_barcodes is not None:
        adata_cs = annotate_doublets(adata_cs, doublet_barcodes)

    # 5. Filter and Save
    generate_visualizations(adata_cs, args.image_prefix)
    compare_umaps(adata_cs, args.image_prefix)

    logger.info(f"Filtering cells (is_empty=False, alpha_hat <= {args.alpha_threshold})")
    keep_mask = (~adata_cs.obs["is_empty"]) & (adata_cs.obs["alpha_hat"] <= args.alpha_threshold)
    adata_filtered = adata_cs[keep_mask].copy()
    generate_noise_boxplot(adata_filtered, args.image_prefix, 'after')

    adata_filtered.write_h5ad(args.cs_filtered_h5ad)
    logger.info("CellSweep complete. Output files written.")


if __name__ == "__main__":
    main()
