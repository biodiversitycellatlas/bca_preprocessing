#!/usr/bin/env python3
"""
Reading the consensus doublet calls that ``COMBINE_DOUBLET_RESULTS`` produces.

Demuxafy's ``Combine_Results.R`` writes one row per evaluated barcode and one
``<method>_DropletType`` column per consensus method; ``params.doublet_consensus_method``
picks which one counts (``AnySinglet`` = doublet only where every tool agrees,
``AnyDoublet`` = doublet where any tool calls it).

Shared by ``filter_doublets.py``, ``run_cellsweep.py`` and ``mtx_to_h5ad.py``, which
all have to agree on that column and on which barcodes it names.
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger("DoubletCalls")

# Each tool's prefix in the Combine_Results.R table, and the obs columns its call and
# score land in
DOUBLET_COL = "doublet_status"
CALLERS = {"scrublet": "scrublet", "scDblFinder": "scdblfinder"}
METHOD_UNS = "doublet_method"


def read_combined_results(combined_results):
    """The Combine_Results.R table, as written."""
    return pd.read_csv(combined_results, sep="\t")


def doublet_barcodes(combined, method, source=""):
    """
    The barcodes `method` calls doublets, out of an already-read table.

    Raises when the method's column is absent: a missing consensus column would
    otherwise read as "no doublets" and pass silently.
    """
    classification_col = f"{method}_DropletType"
    if classification_col not in combined.columns:
        raise ValueError(
            f"Expected column '{classification_col}' not found in {source or 'combined results'}. "
            f"Available columns: {list(combined.columns)}"
        )

    return set(combined.loc[combined[classification_col].str.lower() == "doublet", "Barcode"])


def read_doublet_barcodes(combined_results, method):
    """The barcodes `method` calls doublets, read straight from the file."""
    combined = read_combined_results(combined_results)
    return doublet_barcodes(combined, method, source=combined_results)


def annotate(adata, barcodes):
    """
    Flag the consensus doublets in ``.obs[DOUBLET_COL]``.

    Doublets are called on the cell-called matrix, so only a subset of a raw
    matrix's barcodes was ever evaluated; everything else is left as 'singlet'.
    """
    is_doublet = adata.obs_names.isin(barcodes)
    adata.obs[DOUBLET_COL] = pd.Categorical(
        np.where(is_doublet, "doublet", "singlet"),
        categories=["singlet", "doublet"],
    )
    logger.info(f"Annotated {int(is_doublet.sum())} / {adata.n_obs} barcodes as consensus doublets")
    return adata


def _score_column(combined, prefix):
    """The tool's score column (Scrublet writes ``_Scores``, scDblFinder ``_Score``), or None."""
    return next((c for c in combined.columns
                 if c.lower().startswith(prefix.lower() + "_") and "score" in c.lower()), None)


def annotate_callers(adata, combined, method):
    """
    Each tool's own call and score, beside the consensus: ``.obs['<tool>_call']``
    (singlet / doublet) and ``.obs['<tool>_score']``, NaN for the barcodes never
    evaluated, and ``.uns['doublet_method']``, the consensus rule. A tool missing
    from the table leaves its columns out; a missing score only the score.
    """
    table = combined.set_index("Barcode")
    for prefix, name in CALLERS.items():
        call_col = f"{prefix}_DropletType"
        if call_col not in table.columns:
            logger.warning(f"No '{call_col}' column in the combined results; {name} calls left out")
            continue
        calls = table[call_col].astype(str).str.lower().reindex(adata.obs_names)
        adata.obs[f"{name}_call"] = pd.Categorical(calls, categories=["singlet", "doublet"])
        score_col = _score_column(table, prefix)
        if score_col:
            adata.obs[f"{name}_score"] = pd.to_numeric(table[score_col], errors="coerce").reindex(adata.obs_names).to_numpy()
    adata.uns[METHOD_UNS] = method
    return adata
