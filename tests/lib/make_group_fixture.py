#!/usr/bin/env python3
"""
Two samples of one 'filtering_report' group, for tests/checks/sample_groups.sh.

Sample S1 is make_metacell_fixture.py's sample as it stands. Sample S2 holds the same
cells under other barcodes, except the first SHARED, which both captures happen to
use. Its CellSweep results (alpha_hat, ambient_hat), read metrics and CellReads.stats
are shifted, so the tests can tell which sample a pooled number came from. A third
object, counted against a renamed gene, is the group that must be refused.

Usage: make_group_fixture.py OUTDIR
"""

import os
import sys

import numpy as np
import pandas as pd

import make_metacell_fixture as F

SHARED = 2              # barcodes S2 shares with S1
ALPHA_SHIFT = 0.1       # S2's alpha_hat = S1's + this
AMBIENT_SCALE = 2.0     # S2's ambient_hat = S1's x this
READS_SHIFT = 1000      # S2's counted reads per cell = S1's + this
S2_FRACTIONS = {"mt": 2 * F.FC_MT_FRACTION, "rrna": 2 * F.FC_RRNA_FRACTION, "antisense": 2 * F.FC_ANTISENSE_FRACTION}


def s2_barcodes():
    return [bc if i < SHARED else bc.replace("BC", "XB", 1) for i, bc in enumerate(F.barcodes())]


def s2_reads():
    return [r + READS_SHIFT for r in F.cell_reads()]


def _rename(adata, mapping, sample_id):
    adata = adata.copy()
    adata.obs_names = [mapping.get(bc, bc) for bc in adata.obs_names]
    if "sample_id" in adata.obs.columns:
        adata.obs["sample_id"] = sample_id
    return adata


def write_s2_metrics(outdir):
    with open(os.path.join(outdir, "S2_starsolo_mt_rrna_metrics.txt"), "w") as fh:
        fh.write("Metric,Count\n")
        fh.write(f'"Percentage of rRNA reads (of mapped reads, primary alignment, called cells)",{S2_FRACTIONS["rrna"]:.4f}\n')
        fh.write(f'"Percentage of mtDNA reads (of mapped reads, primary alignment, called cells)",{S2_FRACTIONS["mt"]:.4f}\n')
    with open(os.path.join(outdir, "S2_starsolo_antisense_metrics.txt"), "w") as fh:
        fh.write("Metric,Count\n")
        fh.write(f'Percentage of antisense reads (of reads assigned to genes),{S2_FRACTIONS["antisense"]:.4f}\n')


def write_s2_cellreads(path):
    X = F.counts()
    cols = ["CB", "cbMatch", "cbPerfect", "exonic", "intronic", "mito", "genomicNoGene",
            "countedU", "countedM", "nGenesUnique", "nGenesMulti", "nUMIunique", "nUMImulti"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, (bc, reads) in enumerate(zip(s2_barcodes(), s2_reads())):
            umis, genes = int(X[i].sum()), int((X[i] > 0).sum())
            fh.write(f"{bc}\t{reads}\t{reads}\t{reads}\t0\t0\t0\t{reads}\t0\t{genes}\t0\t{umis}\t0\n")


def main(outdir):
    import anndata as ad

    F.main(outdir)
    cells = ad.read_h5ad(os.path.join(outdir, "cells.h5ad"))
    raw = ad.read_h5ad(os.path.join(outdir, "raw.h5ad"))
    cells.obs["sample_id"] = "S1_starsolo"
    cells.write_h5ad(os.path.join(outdir, "cells_S1.h5ad"))

    mapping = dict(zip(F.barcodes(), s2_barcodes()))
    s2 = _rename(cells, mapping, "S2_starsolo")
    s2.write_h5ad(os.path.join(outdir, "cells_S2.h5ad"))

    s2_raw = _rename(raw, mapping, None)
    s2_raw.obs["alpha_hat"] = s2_raw.obs["alpha_hat"].to_numpy(dtype=float) + ALPHA_SHIFT
    s2_raw.var["ambient_hat"] = s2_raw.var["ambient_hat"].to_numpy(dtype=float) * AMBIENT_SCALE
    s2_raw.write_h5ad(os.path.join(outdir, "raw_S2.h5ad"))

    write_s2_metrics(outdir)
    write_s2_cellreads(os.path.join(outdir, "S2_starsolo_CellReads.stats"))

    # Counted against another annotation: one gene renamed
    bad = s2.copy()
    bad.var_names = ["RENAMED" if g == F.GENE_IDS[-1] else g for g in bad.var_names]
    bad.write_h5ad(os.path.join(outdir, "cells_bad.h5ad"))

    print(f"Group fixture written to {outdir}: S1 and S2, {len(F.barcodes())} cells each, {SHARED} shared barcodes")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main(sys.argv[1])
