#!/usr/bin/env python3
"""
Compute per-cell intronic, mitochondrial, rRNA metrics and generate comparison plots.
"""
import os
import argparse
import json
from collections import defaultdict

import scipy.io
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-cell intronic, mitochondrial, rRNA metrics "
            "and generate comparison plots"
        )
    )
    parser.add_argument(
        "-s", "--solo-output", required=True,
        help="Path to STARsolo output directory (e.g. /.../Sample_Solo.out)"
    )
    parser.add_argument(
        "-m", "--mt-contig", required=True, nargs="+", metavar="CONTIG",
        help=(
            "Name(s) of the mitochondrial contig(s), whitespace separated "
            "(e.g. OW052000.1, or 'chrM M MT'). Reads on any of them count "
            "towards the mitochondrial fraction."
        )
    )
    parser.add_argument(
        "-g", "--gtf", required=True,
        help="Path to the reference GTF file"
    )
    parser.add_argument(
        "-b", "--bam", required=True,
        help="Path to the BAM file"
    )
    parser.add_argument(
        "-o", "--outdir", required=True,
        help="Path to the output directory"
    )
    parser.add_argument(
        "-t", "--min-reads", type=int, required=True,
        help=(
            "Minimum total reads to call a cell (STARsolo "
            "cell/non-cell threshold)"
        )
    )
    parser.add_argument(
        "-w", "--whitelist", required=False,
        help="Path to an optional text file containing rRNA gene names/IDs (one per row)"
    )
    return parser.parse_args()


def load_rrna_intervals(gtf_path, whitelist_path=None):
    """
    Parse GTF file and return a mapping of chromosome -> rRNA intervals.
    """
    # Load whitelist if provided
    whitelist = set()
    if whitelist_path and os.path.exists(whitelist_path):
        with open(whitelist_path) as f:
            whitelist = {line.strip() for line in f if line.strip()}

    intervals = defaultdict(list)
    with open(gtf_path) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            # Fixed GTF column indexing: end is cols[4], attributes are cols[8]
            chrom = cols[0]
            ftype = cols[2]
            start = int(cols[3])
            end = int(cols[4])
            attributes = cols[8]

            is_rrna = False

            # Check if feature type is explicitly rRNA
            if ftype.lower() == "rrna":
                is_rrna = True

            # Check if attributes contain rRNA biotype (handles both GENCODE 'gene_type' and Ensembl 'gene_biotype')
            elif 'gene_biotype "rRNA"' in attributes or 'gene_type "rRNA"' in attributes:
                is_rrna = True
            elif 'transcript_biotype "rRNA"' in attributes or 'transcript_type "rRNA"' in attributes:
                is_rrna = True

            # Check if feature matches a user-provided whitelist
            elif whitelist:
                for item in whitelist:
                    if f'"{item}"' in attributes:
                        is_rrna = True
                        break

            if is_rrna:
                # Convert 1-based inclusive to 0-based half-open
                intervals[chrom].append((start - 1, end))

    return intervals


def load_matrices(solo_dir):
    """
    Load GeneFull_Ex50pAS and Gene matrices and computes intronic percentages.
    Aligns barcodes safely to prevent IndexErrors.
    """
    # File paths
    full_path = os.path.join(solo_dir, "GeneFull_Ex50pAS/raw/matrix.mtx")
    gene_path = os.path.join(solo_dir, "Gene/raw/matrix.mtx")
    full_bc_path = os.path.join(solo_dir, "GeneFull_Ex50pAS/raw/barcodes.tsv")
    gene_bc_path = os.path.join(solo_dir, "Gene/raw/barcodes.tsv")

    # Load matrix data and sum across features (axis=0)
    mat_full = scipy.io.mmread(full_path).tocsc()
    mat_gene = scipy.io.mmread(gene_path).tocsc()

    counts_full_raw = np.array(mat_full.sum(axis=0)).ravel()
    counts_gene_raw = np.array(mat_gene.sum(axis=0)).ravel()

    # Read barcodes using your existing helper function
    bc_full = read_barcodes(full_bc_path)
    bc_gene = read_barcodes(gene_bc_path)

    # Use Pandas to align gene counts to full counts by barcode string
    s_full = pd.Series(counts_full_raw, index=bc_full)
    s_gene = pd.Series(counts_gene_raw, index=bc_gene)

    # Reindex s_gene to match the barcodes of s_full, filling missing ones with 0
    s_gene_aligned = s_gene.reindex(s_full.index, fill_value=0)

    # Extract back to NumPy arrays
    counts_full = s_full.values
    counts_gene = s_gene_aligned.values

    # Avoid division by zero and compute intronic percent
    intronic = np.zeros_like(counts_full, dtype=float)
    mask = counts_full > 0
    intronic[mask] = (
        counts_full[mask] - counts_gene[mask]
    ) / counts_full[mask]

    intronic[intronic < 0] = 0
    return intronic


def read_barcodes(barcodes_path):
    """
    Read barcodes TSV and return a list of cell barcodes.
    """
    df = pd.read_csv(barcodes_path, sep="\t", header=None)
    if df.shape[1] > 1:
        return df.iloc[:, 0].tolist()
    return df.squeeze().tolist()


def scan_bam(
    bam_path, cb_list, mt_contigs, rrna_intervals
):
    """
    Scan BAM file to count total, mitochondrial, and rRNA reads per barcode.
    A read counts as mitochondrial when it maps to any of mt_contigs.
    Returns three dicts: total_counts, mt_counts, rrna_counts.
    """
    total = defaultdict(int)
    mt = defaultdict(int)
    rrna = defaultdict(int)

    mt_set = set(mt_contigs)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # The default contig names are a guess covering common assemblies, so
        # report which ones actually exist rather than silently reporting 0%
        present = mt_set.intersection(bam.references)
        if present:
            print(f"Mitochondrial contigs found in BAM: {' '.join(sorted(present))}")
        else:
            print(
                f"WARNING: none of the requested mitochondrial contigs "
                f"({' '.join(sorted(mt_set))}) are present in {bam_path}; "
                "mtDNA percentages will be 0"
            )

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            try:
                cb = read.get_tag("CB")
            except KeyError:
                continue
            total[cb] += 1
            if read.reference_name in mt_set:
                mt[cb] += 1
            chrom = read.reference_name
            if chrom in rrna_intervals:
                rstart, rend = read.reference_end, read.reference_start
                for istart, iend in rrna_intervals[chrom]:
                    # Overlap check
                    if istart < rstart and iend > rend:
                        rrna[cb] += 1
                        break
                    # Broad check for any overlap
                    if max(rstart, istart) < min(rend, iend):
                         rrna[cb] += 1
                         break

    # Build arrays aligned to barcode list
    total_arr = np.array([total.get(c, 0) for c in cb_list])
    mt_arr = np.array([mt.get(c, 0) for c in cb_list])
    rrna_arr = np.array([rrna.get(c, 0) for c in cb_list])
    return total_arr, mt_arr, rrna_arr


def compute_percentages(total, part):
    """
    Compute fraction part/total, handling zeros.
    """
    pct = np.zeros_like(part, dtype=float)
    mask = total > 0
    pct[mask] = part[mask] / total[mask]
    return pct


def save_metrics(df, out_dir, prefix):
    """
    Save DataFrame as CSV and JSON for interactive dashboards.
    JSON is saved in column-oriented format (dict of lists) for efficiency.
    """
    csv_path = os.path.join(out_dir, f"{prefix}_metrics.csv")
    json_path = os.path.join(out_dir, f"{prefix}_metrics.json")

    # Save CSV
    df.to_csv(csv_path, index=False)

    # Filter out zero-read cells for JSON
    df_json = df[df["TotalReads"] > 0].copy()

    # Round floats to 4 decimals to save space
    float_cols = ["IntronicPercent", "MTPercent", "rRNAPercent"]
    df_json[float_cols] = df_json[float_cols].round(4)

    # Convert to dictionary of lists {"col": [val, val...]}
    metrics_dict = df_json.to_dict(orient="list")

    # Save using json format
    with open(json_path, "w") as f:
        json.dump(metrics_dict, f)

    print(f"Saved metrics CSV to {csv_path}")
    print(f"Saved optimized JSON to {json_path}")


def plot_comparisons(
    df, out_dir, prefix, min_reads
):
    """
    Generate and save scatter plots for predefined metric pairs.
    """
    axis_labels = {
        "IntronicPercent": "Intronic Reads (%)",
        "MTPercent": "mtDNA (%)",
        "rRNAPercent": "rRNA (%)",
        "TotalReads": "Total Reads (Cell Size)",
    }
    comparisons = [
        ("IntronicPercent", "MTPercent", "Percentage Intronic vs Mitochondrial Reads Per Cell"),
        ("IntronicPercent", "rRNAPercent", "Percentage Intronic vs rRNA Reads Per Cell"),
        ("IntronicPercent", "TotalReads", "Percentage Intronic Reads Per Cell vs Cell Size"),
        ("MTPercent", "rRNAPercent", "Percentage Mitochondrial vs rRNA Reads Per Cell"),
        ("MTPercent", "TotalReads", "Percentage Mitochondrial Reads Per Cell vs Cell Size"),
        ("rRNAPercent", "TotalReads", "Percentage rRNA Reads Per Cell vs Cell Size"),
    ]
    low_mask = df["TotalReads"] < min_reads

    for x, y, title in comparisons:
        fig, ax = plt.subplots()
        ax.scatter(
            df.loc[low_mask, x],
            df.loc[low_mask, y],
            s=5, color="#b6b5b5",
            label=f"< {min_reads} reads"
        )
        ax.scatter(
            df.loc[~low_mask, x],
            df.loc[~low_mask, y],
            s=5, color="steelblue",
            label=f"≥ {min_reads} reads"
        )
        ax.set_xlabel(axis_labels[x])
        ax.set_ylabel(axis_labels[y])
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

    # Load data
    intronic = load_matrices(solo_out)
    barcodes = read_barcodes(
        os.path.join(
            solo_out, "GeneFull_Ex50pAS/raw/barcodes.tsv"
        )
    )
    rrna_int = load_rrna_intervals(args.gtf, args.whitelist)

    # Compute read counts
    total_cnt, mt_cnt, rrna_cnt = scan_bam(
        args.bam, barcodes, args.mt_contig, rrna_int
    )

    # Calculate percentages
    mt_pct = compute_percentages(total_cnt, mt_cnt)
    rrna_pct = compute_percentages(total_cnt, rrna_cnt)

    # Build DataFrame
    df = pd.DataFrame({
        "Cell": barcodes,
        "IntronicPercent": intronic,
        "MTPercent": mt_pct,
        "rRNAPercent": rrna_pct,
        "TotalReads": total_cnt,
    })

    # Save metrics and plots
    save_metrics(df, args.outdir, prefix)
    plot_comparisons(df, args.outdir, prefix, args.min_reads)


if __name__ == "__main__":
    main()
