#!/usr/bin/env python3
"""
Build a synthetic sample whose metacell assignment is planted, for
tests/checks/metacell_filtering.sh.

Metacell2 itself is not needed: the cells h5ad carries the planted assignment in
``obs['mc_truth']``, which ``run_metacells.py --metacell_obs_col mc_truth`` takes
instead of grouping. Everything downstream of the grouping -- the per-metacell QC,
the gene matching, the CellSweep lift, the report and the apply step -- then runs
on numbers known by construction.

    group      cells   per-cell counts                    doublets    alpha_hat
    M0.1       30      PROFILES['M0.1'] x (1 + c % 3)     cells 0-2   0.01
    M1.2       30      mito-rich                          none        0.02
    M2.3       30      ...                                cells 0-19  0.03
    __outliers__  5
    __excluded__  5

Each cell's counts are its group's profile times (1 + c % 3), c the cell's index
within the group, so the per-cell mito percentages vary within a group and the
pooled one is not simply their mean. Expected values are recomputed from the dense
matrix, never through the code under test.

The gene axis carries the awkward cases the matching has to handle:

    MT-CO1, MT-ND1          on chrM: mitochondrial
    RRNA1                   biotype rRNA
    GENE003                 domains keyed by transcript id, with a version suffix
    GENE004                 domains keyed by protein id, accession with a version
    LOC5-AS1                a gene id with '-', which alevin-fry reduces to 'LOC5'
    GENE009                 no UMI anywhere
    GENE010                 named '</script><script>alert(1)'
    GENE011                 one UMI, named like the report's placeholder

The ambient (raw) h5ad holds the called cells shuffled among empty droplets, so a
lift by position instead of by barcode gives wrong values.

Usage: make_metacell_fixture.py OUTDIR
"""

import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

# gene_id, gene_name, chrom, biotype
GENES = [
    ("GENE001", "Actb", "chr1", "protein_coding"),
    ("GENE002", "Gapdh", "chr1", "protein_coding"),
    ("GENE003", "Rpl3", "chr2", "protein_coding"),
    ("GENE004", "H2a", "chr2", "protein_coding"),
    ("LOC5-AS1", "", "chr3", "lncRNA"),
    ("MT-CO1", "mt-Co1", "chrM", "protein_coding"),
    ("MT-ND1", "mt-Nd1", "chrM", "protein_coding"),
    ("RRNA1", "Rn45s", "chr4", "rRNA"),
    ("GENE009", "Silent", "chr5", "protein_coding"),
    ("GENE010", "</script><script>alert(1)", "chr5", "protein_coding"),
    ("GENE011", "__REPORT_DATA_PLACEHOLDER__", "chr5", "protein_coding"),
]
GENE_IDS = [g[0] for g in GENES]
MITO_GENES = ["MT-CO1", "MT-ND1"]
RRNA_GENES = ["RRNA1"]
ZERO_GENE = "GENE009"
SCRIPT_GENE = "GENE010"
LOW_GENE = "GENE011"

# Expected domains after build_gene_table.py has resolved and cleaned them
EXPECTED_DOMAINS = {"GENE003": ["Ribosomal_L3"], "GENE004": ["PF00125"]}

MT_CONTIG = ["chrM", "MT"]

GROUPS = ["M0.1", "M1.2", "M2.3"]
CELLS_PER_GROUP = 30
N_OUTLIERS = 5
N_EXCLUDED = 5
N_EMPTY = 20

# Counts per gene, in GENE_IDS order (GENE009-011 are set separately)
PROFILES = {
    "M0.1":         [20, 10, 8, 2, 3, 2, 1, 1],
    "M1.2":         [5, 20, 4, 6, 0, 10, 5, 0],
    "M2.3":         [10, 10, 10, 0, 5, 1, 0, 4],
    "__outliers__": [3, 3, 0, 0, 0, 0, 0, 0],
    "__excluded__": [1, 1, 0, 0, 0, 0, 0, 0],
}
DOUBLETS = {"M0.1": 3, "M1.2": 0, "M2.3": 20, "__outliers__": 0, "__excluded__": 0}
ALPHA = {"M0.1": 0.01, "M1.2": 0.02, "M2.3": 0.03, "__outliers__": 0.5, "__excluded__": 0.9}

# The report shows the called-cell values; the whole-BAM ones differ, so showing
# them instead would fail the assertion
FC_MT_FRACTION = 0.0123
FC_RRNA_FRACTION = 0.0456
FC_MT_BAM_FRACTION = 0.0321
FC_RRNA_BAM_FRACTION = 0.0654
FC_ANTISENSE_FRACTION = 0.0789


def group_sizes():
    sizes = {g: CELLS_PER_GROUP for g in GROUPS}
    sizes["__outliers__"] = N_OUTLIERS
    sizes["__excluded__"] = N_EXCLUDED
    return sizes


def cells():
    """(barcode, group, index within group), in the order the cells h5ad holds them."""
    out = []
    n = 0
    for g, size in group_sizes().items():
        for c in range(size):
            out.append((f"BC{n:04d}", g, c))
            n += 1
    return out


def barcodes():
    return [bc for bc, _, _ in cells()]


def groups():
    return [g for _, g, _ in cells()]


def counts():
    """Cells x genes raw counts, dense."""
    rows = []
    for _, g, c in cells():
        row = [v * (1 + c % 3) for v in PROFILES[g]] + [0, 0, 0]
        row[GENE_IDS.index(SCRIPT_GENE)] = 1 if c % 5 == 0 else 0
        rows.append(row)
    X = np.array(rows, dtype=np.float64)
    X[0, GENE_IDS.index(LOW_GENE)] = 1
    return X


def doublet_status():
    return ["doublet" if c < DOUBLETS[g] else "singlet" for _, g, c in cells()]


def alpha_hat():
    return [ALPHA[g] for g in groups()]


def ambient_hat(j):
    return round(0.001 * (j + 1), 6)


def denoised(X):
    """CellSweep's denoised counts: deterministic, and different from the raw ones."""
    return np.floor(X / 2)


def _mito_columns():
    return [GENE_IDS.index(g) for g in MITO_GENES]


def expected_mito_pct(group):
    """Pooled over the group's cells, from the dense counts: mito UMIs over all UMIs."""
    X = counts()
    rows = [i for i, g in enumerate(groups()) if g == group]
    total = X[rows].sum()
    return 100.0 * X[np.ix_(rows, _mito_columns())].sum() / total if total else 0.0


def expected_cell_mito_pct():
    X = counts()
    total = X.sum(axis=1)
    return np.where(total > 0, 100.0 * X[:, _mito_columns()].sum(axis=1) / np.where(total > 0, total, 1), 0.0)


def expected_doublet_pct(group):
    return 100.0 * DOUBLETS[group] / group_sizes()[group]


def write_gtf(path):
    with open(path, "w") as fh:
        fh.write("#!genome-build fixture\n")
        for i, (gid, name, chrom, biotype) in enumerate(GENES):
            start, end = 1000 * (i + 1), 1000 * (i + 1) + 500
            name_attr = f' gene_name "{name}";' if name else ""
            fh.write(f'{chrom}\tfixture\tgene\t{start}\t{end}\t.\t+\t.\tgene_id "{gid}";{name_attr} gene_biotype "{biotype}";\n')
            tx = f"TX{i + 1:03d}"
            fh.write(f'{chrom}\tfixture\ttranscript\t{start}\t{end}\t.\t+\t.\tgene_id "{gid}"; transcript_id "{tx}";{name_attr} gene_biotype "{biotype}";\n')
            fh.write(f'{chrom}\tfixture\tCDS\t{start}\t{end}\t.\t+\t0\tgene_id "{gid}"; transcript_id "{tx}"; protein_id "PROT{i + 1:03d}";\n')


def write_annotation(path):
    # eggNOG-mapper style: '#query' header, transcript/protein keys, '-' for none
    rows = [
        ("TX003.1", "Ribosomal_L3"),      # transcript id, version suffix on the key
        ("PROT004", "PF00125.27"),        # protein id, versioned accession
        ("GENE001", "-"),                 # resolves, but no domain
        ("UNKNOWN_KEY", "Foo"),           # resolves to nothing
    ]
    with open(path, "w") as fh:
        fh.write("#query\tPFAMs\n")
        for key, dom in rows:
            fh.write(f"{key}\t{dom}\n")


def write_metrics(outdir, sample_id):
    with open(os.path.join(outdir, f"{sample_id}_mt_rrna_metrics.txt"), "w") as fh:
        fh.write("Metric,Count\n")
        # Metric names carry commas and are quoted, as calculate_read_metrics.py writes them
        fh.write(f'"Percentage of rRNA reads (of mapped reads, primary alignment)",{FC_RRNA_BAM_FRACTION:.4f}\n')
        fh.write(f'"Percentage of mtDNA reads (of mapped reads, primary alignment)",{FC_MT_BAM_FRACTION:.4f}\n')
        fh.write(f'"Percentage of rRNA reads (of mapped reads, primary alignment, called cells)",{FC_RRNA_FRACTION:.4f}\n')
        fh.write(f'"Percentage of mtDNA reads (of mapped reads, primary alignment, called cells)",{FC_MT_FRACTION:.4f}\n')
    with open(os.path.join(outdir, f"{sample_id}_antisense_metrics.txt"), "w") as fh:
        fh.write("Metric,Count\n")
        fh.write(f"Percentage of antisense reads (of reads assigned to genes),{FC_ANTISENSE_FRACTION:.4f}\n")


SMOKE_TYPES = 4
SMOKE_CELLS_PER_TYPE = 150
SMOKE_GENES = 300


def write_mc2_smoke(outdir):
    """
    A dataset Metacell2 can actually group, for the mc2_smoke case: four cell types,
    each with its own 30 marker genes over a shared background, Poisson counts. The
    planted fixture above is far too small for Metacell2 to select any gene on.
    """
    import anndata as ad

    rng = np.random.default_rng(11)
    # ~5,000 UMIs per cell: Metacell2 also sizes metacells by UMIs, so shallower cells
    # give too few metacells to lay out
    base = 8 * rng.gamma(0.6, 1.0, SMOKE_GENES)
    rows, labels = [], []
    for t in range(SMOKE_TYPES):
        profile = base.copy()
        profile[20 + 30 * t: 50 + 30 * t] *= 25
        profile[:2] = 24.0                                  # two mitochondrial genes
        for _ in range(SMOKE_CELLS_PER_TYPE):
            rows.append(rng.poisson(profile * rng.uniform(0.7, 1.4)))
            labels.append(f"type{t}")
    ids = [f"MT-S{j}" if j < 2 else f"SMK{j:04d}" for j in range(SMOKE_GENES)]
    obs = pd.DataFrame({"truth": labels}, index=pd.Index([f"SC{i:05d}" for i in range(len(rows))]))
    ad.AnnData(X=sp.csr_matrix(np.array(rows, dtype=np.float32)), obs=obs,
               var=pd.DataFrame(index=pd.Index(ids))).write_h5ad(os.path.join(outdir, "mc2_cells.h5ad"))

    table = pd.DataFrame({
        "gene_id": ids, "gene_name": ids, "chrom": ["chrM" if j < 2 else "chr1" for j in range(SMOKE_GENES)],
        "biotype": "protein_coding", "is_mito": [int(j < 2) for j in range(SMOKE_GENES)], "is_rrna": 0, "pfam": "",
    })
    table.to_csv(os.path.join(outdir, "mc2_gene_table.tsv"), sep="\t", index=False)


def main(outdir):
    import anndata as ad

    os.makedirs(outdir, exist_ok=True)
    write_gtf(os.path.join(outdir, "genes.gtf"))
    write_annotation(os.path.join(outdir, "pfam.tsv"))
    write_metrics(outdir, "S1_starsolo")

    X = counts()
    bcs = barcodes()
    var = pd.DataFrame(index=pd.Index(GENE_IDS))

    # The cell-called object, as MTX_TO_H5AD writes it, plus the planted assignment
    obs = pd.DataFrame(index=pd.Index(bcs))
    obs["sample_id"] = "S1"
    obs["doublet_status"] = pd.Categorical(doublet_status(), categories=["singlet", "doublet"])
    obs["mc_truth"] = groups()
    ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var.copy()).write_h5ad(os.path.join(outdir, "cells.h5ad"))

    # The same without doublet calls, for the 'absent' mode
    ad.AnnData(X=sp.csr_matrix(X), obs=obs.drop(columns=["doublet_status"]), var=var.copy()).write_h5ad(
        os.path.join(outdir, "cells_nodoublets.h5ad"))

    # The raw object: called cells shuffled among empty droplets, CellSweep's results on it
    rng = np.random.default_rng(7)
    empty = [f"EMPTY{i:03d}" for i in range(N_EMPTY)]
    raw_bcs = bcs + empty
    raw_X = np.vstack([X, np.zeros((N_EMPTY, X.shape[1]))])
    raw_X[len(bcs):, 0] = 1
    order = rng.permutation(len(raw_bcs))
    raw_bcs = [raw_bcs[i] for i in order]
    raw_X = raw_X[order]

    raw_obs = pd.DataFrame(index=pd.Index(raw_bcs))
    lookup = dict(zip(bcs, alpha_hat()))
    raw_obs["alpha_hat"] = [lookup.get(bc, 0.99) for bc in raw_bcs]
    raw_obs["is_empty"] = [bc.startswith("EMPTY") for bc in raw_bcs]
    raw_obs["celltype"] = "fixture"
    raw_var = var.copy()
    raw_var["ambient_hat"] = [ambient_hat(j) for j in range(len(GENE_IDS))]
    raw = ad.AnnData(X=sp.csr_matrix(raw_X), obs=raw_obs, var=raw_var)
    raw.layers["cellsweep"] = sp.csr_matrix(denoised(raw_X))
    raw.write_h5ad(os.path.join(outdir, "raw.h5ad"))

    write_mc2_smoke(outdir)

    print(f"Fixture written to {outdir}: {len(bcs)} called cells, {len(GENE_IDS)} genes")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
