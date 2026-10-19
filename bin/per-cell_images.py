#!/usr/bin/env python3
"""
Compute per-cell intronic, mtDNA and rRNA metrics and generate comparison plots.

The read counts come from ``<id>_barcode_reads.tsv.gz``, written by
calculate_read_metrics.py (CALC_READ_METRICS) in its pass over the BAM, so this
script never reads the BAM itself. Every read metric counts a read once, at its
primary alignment, and matches the called-cell rows of ``_mt_rrna_metrics.txt``:

    MappedReads       reads with a CB tag (primary alignment)
    MTPercent         % of MappedReads whose primary alignment is on a mitochondrial contig
    rRNAPercent       % of MappedReads overlapping an rRNA region (any rRNA biotype, plus
                      the whole of every contig of the added rRNA reference)
    IntronicPercent   % of uniquely mapped reads STARsolo flags intronic (CellReads.stats)
    UnsplicedPercent  % of Velocyto UMIs classed unspliced, when Velocyto was run

All percentages are written as 0-100.
"""
import os
import argparse
import json
import sys

import scipy.io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-cell intronic, mtDNA and rRNA metrics "
            "and generate comparison plots"
        )
    )
    parser.add_argument(
        "-s", "--solo-output", required=True,
        help="Path to STARsolo output directory (e.g. /.../Sample_Solo.out)"
    )
    parser.add_argument(
        "-r", "--read-counts", required=True,
        help=(
            "Mapped, mtDNA and rRNA reads per barcode (<id>_barcode_reads.tsv.gz "
            "from calculate_read_metrics.py)"
        )
    )
    parser.add_argument(
        "-o", "--outdir", required=True,
        help="Path to the output directory"
    )
    parser.add_argument(
        "-t", "--min-reads", type=int, default=0,
        help=(
            "Minimum mapped reads to call a cell. Only used as a fallback when "
            "--cell-barcodes is not given"
        )
    )
    parser.add_argument(
        "-c", "--cell-barcodes", required=False,
        help=(
            "Path to the barcodes.tsv of the filtered count matrix. When given, "
            "cells are split on membership of that set instead of on --min-reads, "
            "which makes the plots agree with the matrix that was actually "
            "filtered"
        )
    )
    return parser.parse_args()


def load_cell_barcodes(path):
    """Read the called-cell barcodes from a filtered matrix's ``barcodes.tsv``.

    Returns a set of barcodes, or ``None`` when no usable file was given, in
    which case the caller falls back to the read threshold.
    """
    if not path or not os.path.exists(path):
        return None

    try:
        df = pd.read_csv(path, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        # A cell filter that kept nothing would have failed upstream, so an empty
        # file means the list is unusable rather than that there are no cells
        print(f"Warning: {path} is empty; falling back to the read threshold")
        return None

    if df.empty:
        return None

    return set(df.iloc[:, 0].astype(str))


def load_cellreads_intronic(solo_dir):
    """
    Per-barcode intronic read percentage from STARsolo's ``CellReads.stats``.

    ``intronic / genomeU``: STARsolo flags reads intronic only when they have a
    unique genomic alignment, so ``genomeU`` is the denominator, the same one the
    dashboard and mapping_stats.tsv use for their library- and cell-level values.
    Returns a Series indexed by barcode, or ``None`` when the file is absent.
    """
    path = os.path.join(solo_dir, "GeneFull_Ex50pAS", "CellReads.stats")
    if not os.path.exists(path):
        print(
            f"Warning: no {path}; the intronic percentage needs soloCellReadStats "
            "and is skipped for this sample.",
            file=sys.stderr,
        )
        return None

    stats = pd.read_csv(path, sep="\t", index_col=0)
    if "intronic" not in stats.columns or "genomeU" not in stats.columns:
        print(f"Warning: 'intronic'/'genomeU' columns not found in {path}", file=sys.stderr)
        return None
    # Aggregate row for reads whose barcode is not in the passlist
    stats = stats.drop(index="CBnotInPasslist", errors="ignore")
    stats.index = stats.index.astype(str)
    return compute_percentages(stats["genomeU"].to_numpy(), stats["intronic"].to_numpy(),
                               index=stats.index)


def load_velocyto_unspliced(solo_dir):
    """
    Per-barcode unspliced UMI percentage from STARsolo's Velocyto matrices,
    ``unspliced / (spliced + unspliced + ambiguous)``.

    Velocyto classifies each UMI once against the transcript models, so this is the
    UMI-level counterpart of the intronic read percentage. It only covers the Gene
    feature (annotated genes, sense strand), so the two are not expected to be
    equal. Returns a Series indexed by barcode, or ``None`` when Velocyto was not run.
    """
    vdir = os.path.join(solo_dir, "Velocyto", "raw")
    paths = {kind: os.path.join(vdir, f"{kind}.mtx") for kind in ("spliced", "unspliced", "ambiguous")}
    bc_path = os.path.join(vdir, "barcodes.tsv")
    if not all(os.path.exists(p) for p in [*paths.values(), bc_path]):
        return None

    totals = {kind: np.asarray(scipy.io.mmread(p).tocsc().sum(axis=0)).ravel()
              for kind, p in paths.items()}
    barcodes = pd.Index(read_barcodes(bc_path)).astype(str)
    every = totals["spliced"] + totals["unspliced"] + totals["ambiguous"]
    return compute_percentages(every, totals["unspliced"], index=barcodes)


def read_barcodes(barcodes_path):
    """
    Read barcodes TSV and return a list of cell barcodes.
    """
    df = pd.read_csv(barcodes_path, sep="\t", header=None)
    if df.shape[1] > 1:
        return df.iloc[:, 0].tolist()
    return df.squeeze().tolist()


def load_read_counts(path, cb_list):
    """
    Mapped, mtDNA and rRNA reads per barcode from CALC_READ_METRICS'
    ``<id>_barcode_reads.tsv.gz``, as three arrays aligned to *cb_list*.

    Barcodes absent from the file had no read at a primary alignment and get 0;
    barcodes outside *cb_list* (not in the raw matrix) are dropped.
    """
    # keep_default_na: a barcode such as "NA" or "-" is a barcode, not a missing value
    counts = pd.read_csv(path, sep="\t", index_col="CB", dtype={"CB": str}, keep_default_na=False)
    counts = counts.reindex(cb_list, fill_value=0)
    return (counts["MappedReads"].to_numpy(),
            counts["MTReads"].to_numpy(),
            counts["rRNAReads"].to_numpy())


def compute_percentages(total, part, index=None):
    """
    Compute 100 * part / total, with 0 where total is 0. Returns a Series when
    *index* is given, otherwise an array.
    """
    total = np.asarray(total, dtype=float)
    part = np.asarray(part, dtype=float)
    pct = np.zeros_like(part, dtype=float)
    mask = total > 0
    pct[mask] = 100.0 * part[mask] / total[mask]
    return pd.Series(pct, index=index) if index is not None else pct


def save_metrics(df, out_dir, prefix):
    """
    Save DataFrame as CSV and JSON for interactive dashboards.
    JSON is saved in column-oriented format (dict of lists) for efficiency.

    ``IsCell`` travels with the metrics so the dashboard splits its interactive
    plots on the same cell set as the static PNGs.
    """
    csv_path = os.path.join(out_dir, f"{prefix}_metrics.csv")
    json_path = os.path.join(out_dir, f"{prefix}_metrics.json")

    # Save CSV
    df.to_csv(csv_path, index=False)

    # Filter out zero-read cells for JSON
    df_json = df[df["MappedReads"] > 0].copy()

    # Round floats to save space. IntronicPercent and UnsplicedPercent are absent
    # when their inputs were not produced, so the columns are intersected with the frame.
    float_cols = [c for c in PCT_COLUMNS if c in df_json.columns]
    df_json[float_cols] = df_json[float_cols].round(2)

    # Convert to dictionary of lists {"col": [val, val...]}
    metrics_dict = df_json.to_dict(orient="list")

    # Save using json format
    with open(json_path, "w") as f:
        json.dump(metrics_dict, f)

    print(f"Saved metrics CSV to {csv_path}")
    print(f"Saved optimized JSON to {json_path}")


PCT_COLUMNS = ("IntronicPercent", "UnsplicedPercent", "MTPercent", "rRNAPercent")

AXIS_LABELS = {
    "IntronicPercent": "Intronic reads (% of uniquely mapped reads)",
    "UnsplicedPercent": "Unspliced UMIs (%, velocyto)",
    "MTPercent": "mtDNA reads (% of mapped reads)",
    "rRNAPercent": "rRNA reads (% of mapped reads)",
    "MappedReads": "Mapped reads",
}


def plot_comparisons(
    df, out_dir, prefix, cell_label, noncell_label
):
    """
    Generate and save scatter plots for predefined metric pairs, split on the
    ``IsCell`` column.
    """
    comparisons = [
        ("IntronicPercent", "MTPercent", "% Intronic vs % mtDNA Reads Per Cell"),
        ("IntronicPercent", "rRNAPercent", "% Intronic vs % rRNA Reads Per Cell"),
        ("IntronicPercent", "MappedReads", "% Intronic Reads vs Mapped Reads Per Cell"),
        ("MTPercent", "rRNAPercent", "% mtDNA vs % rRNA Reads Per Cell"),
        ("MTPercent", "MappedReads", "% mtDNA Reads vs Mapped Reads Per Cell"),
        ("rRNAPercent", "MappedReads", "% rRNA Reads vs Mapped Reads Per Cell"),
        ("IntronicPercent", "UnsplicedPercent", "% Intronic Reads vs % Unspliced UMIs Per Cell"),
    ]
    # Pairs involving a metric that could not be computed are dropped rather than
    # plotted empty: without CellReads.stats that removes the intronic panels, and
    # without Velocyto the unspliced one.
    comparisons = [(x, y, title) for x, y, title in comparisons
                   if x in df.columns and y in df.columns]

    low_mask = ~df["IsCell"].astype(bool)

    for x, y, title in comparisons:
        fig, ax = plt.subplots()
        ax.scatter(
            df.loc[low_mask, x],
            df.loc[low_mask, y],
            s=5, color="#b6b5b5",
            label=noncell_label
        )
        ax.scatter(
            df.loc[~low_mask, x],
            df.loc[~low_mask, y],
            s=5, color="steelblue",
            label=cell_label
        )
        ax.set_xlabel(AXIS_LABELS[x])
        ax.set_ylabel(AXIS_LABELS[y])
        ax.set_title(title)
        ax.legend(markerscale=2, fontsize="small", framealpha=0.8)
        fig.tight_layout()
        out_png = os.path.join(
            out_dir, f"{prefix}_{x}_vs_{y}.png"
        )
        fig.savefig(out_png)
        plt.close(fig)
        print(f"Saved plot {out_png}")


def main():
    args = parse_args()
    solo_out = args.solo_output.rstrip("/")
    prefix = os.path.basename(solo_out).replace("_Solo.out", "")

    # Every barcode STARsolo saw, in the raw matrix's order
    barcodes = [str(bc) for bc in read_barcodes(
        os.path.join(solo_out, "GeneFull_Ex50pAS/raw/barcodes.tsv")
    )]

    # Read counts per barcode, from CALC_READ_METRICS' pass over the BAM
    total_cnt, mt_cnt, rrna_cnt = load_read_counts(args.read_counts, barcodes)

    # Calculate percentages
    mt_pct = compute_percentages(total_cnt, mt_cnt)
    rrna_pct = compute_percentages(total_cnt, rrna_cnt)

    # Called cells come from the filtered matrix when available. The threshold
    # fallback compares a UMI cutoff against per-barcode read counts, so it only
    # approximates the cell set; the barcode list is the same one the matrix was
    # filtered on and needs no such assumption.
    cell_barcodes = load_cell_barcodes(args.cell_barcodes)
    if cell_barcodes is not None:
        is_cell = np.array([bc in cell_barcodes for bc in barcodes])
        cell_label = f"Called cell (n={int(is_cell.sum())})"
        noncell_label = "Not called"
        print(
            f"Splitting on {len(cell_barcodes)} barcodes from {args.cell_barcodes}; "
            f"{int(is_cell.sum())} of {len(barcodes)} matched"
        )
    else:
        is_cell = total_cnt >= args.min_reads
        cell_label = f"≥ {args.min_reads} reads"
        noncell_label = f"< {args.min_reads} reads"
        print(f"No cell barcode list given; splitting on {args.min_reads} reads")

    # Build DataFrame. The intronic and unspliced columns are omitted entirely
    # rather than filled with a placeholder when their inputs were not produced, so
    # nothing downstream can mistake a stand-in value for a measurement. Barcodes
    # absent from those inputs had no reads there and get 0, like the other columns.
    columns = {"Cell": barcodes}
    intronic = load_cellreads_intronic(solo_out)
    if intronic is not None:
        columns["IntronicPercent"] = intronic.reindex(barcodes, fill_value=0.0).to_numpy()
    unspliced = load_velocyto_unspliced(solo_out)
    if unspliced is not None:
        columns["UnsplicedPercent"] = unspliced.reindex(barcodes, fill_value=0.0).to_numpy()
    columns.update({
        "MTPercent": mt_pct,
        "rRNAPercent": rrna_pct,
        "MappedReads": total_cnt,
        "IsCell": is_cell.astype(int),
    })
    df = pd.DataFrame(columns)

    # Save metrics and plots
    save_metrics(df, args.outdir, prefix)
    plot_comparisons(df, args.outdir, prefix, cell_label, noncell_label)


if __name__ == "__main__":
    main()
