#!/usr/bin/env python3
"""
Reading the consensus doublet calls that ``COMBINE_DOUBLET_RESULTS`` produces.

Demuxafy's ``Combine_Results.R`` writes one row per evaluated barcode and one
``<method>_DropletType`` column per consensus method; ``params.doublet_consensus_method``
picks which one counts (``AnySinglet`` = doublet only where every tool agrees).

Shared by ``filter_doublets.py``, ``run_cellsweep.py`` and ``mtx_to_h5ad.py``, which
all have to agree on that column and on which barcodes it names.
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger("DoubletCalls")

DOUBLET_COL = "doublet_status"


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
