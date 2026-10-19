#!/usr/bin/env python3
"""
Build synthetic ``UMIperCellSorted.txt`` curves with a known cell/ambient boundary.

Cell calling fails silently: a cutoff is always a plausible-looking integer, and
nothing about a wrong one raises.  The only way to test an estimator is to plant
a boundary and check that it is found, so each archetype below is a mixture of
lognormal populations whose sizes *are* the answer -- the assertions restate
``PLANTED_BOUNDARY`` rather than hard-coding a rank some earlier run happened to
produce.

Three archetypes, each targeting a different failure:

``sharp``
    One cell population an order of magnitude above the ambient cloud.  Every
    estimator should find this; one that does not is broken, and the arms should
    agree with each other to within noise.

``two_knee``
    A bright sub-population above the real cells, so the criterion curve has two
    downward bends.  The bright/normal transition at rank 200 is by far the
    deeper of the two -- measured at some 430 noise scales against 96 for the
    cell/ambient boundary at rank 1000 -- because relative barcode weight is
    larger at low rank.  A rule that takes the deepest bend therefore returns
    rank 200, and the biologically correct answer is 1000.  This is the case the
    prominence-plus-nearest-expected selection exists for.

    The sub-population is deliberately sized so that rank 200 falls *inside* the
    0.1x-10x window around 1000 expected cells.  An earlier version put it at
    rank 61, where the window excluded it and both selection rules agreed for
    the wrong reason -- the fixture looked like it exercised the choice and did
    not.

``no_knee``
    A single lognormal with no ambient step at all.  There is no boundary to
    find, and the honest outcome is a loosely determined one: the criterion's
    minimum is broad, so its basin should be far wider than on ``sharp``.  A
    method that reports a confident cutoff here is lying.

Importable, so the checks name the same constants that generated the data.
"""

import argparse
import os
from typing import Dict, Optional, Tuple

import numpy as np


SEED = 20260904
MIN_UMIS = 100

# (count, median UMIs, sigma of log counts) per population, in descending order
# of brightness. Dict order is the draw order, which is what makes the fixtures
# reproducible.
ARCHETYPES: Dict[str, Dict[str, Tuple[int, float, float]]] = {
    "sharp": {
        "cells": (800, 3000.0, 0.50),
        "ambient": (30000, 60.0, 0.65),
    },
    "two_knee": {
        "bright": (200, 30000.0, 0.35),
        "cells": (800, 1500.0, 0.45),
        "ambient": (30000, 60.0, 0.65),
    },
    "no_knee": {
        "all": (20000, 300.0, 0.80),
    },
}

# The rank the cell/ambient boundary sits at: every non-ambient barcode.
PLANTED_BOUNDARY: Dict[str, Optional[int]] = {
    "sharp": 800,
    "two_knee": 1000,
    "no_knee": None,
}

# The competing, deeper bend in ``two_knee``: where the bright sub-population
# gives way to the ordinary cells. A rule that takes the deepest minimum lands
# here instead of on PLANTED_BOUNDARY.
DECOY_BEND: Dict[str, Optional[int]] = {
    "sharp": None,
    "two_knee": 200,
    "no_knee": None,
}

# What the pipeline would be told to expect for each fixture.
EXPECTED_CELLS: Dict[str, Optional[int]] = {
    "sharp": 800,
    "two_knee": 1000,
    "no_knee": 2000,
}

# Half-width of the acceptance band around a planted boundary, in decades of
# rank. 0.15 decades is a factor of 1.41 -- loose enough that the population
# overlap and the smoothing width do not matter, tight enough to tell rank 1000
# from the competing bend at rank 200 in ``two_knee``.
TOLERANCE_DECADES = 0.15

# A curve too short for any estimator to say anything about.
TINY_COUNTS = [4200, 3100, 900, 400, 150]


def lognormal_counts(rng: np.random.Generator, n: int, median: float, sigma: float) -> np.ndarray:
    """One population of per-barcode counts, floored at 1."""
    draws = np.exp(np.log(median) + sigma * rng.standard_normal(n))
    return np.maximum(np.round(draws), 1).astype(np.int64)


def build(kind: str) -> np.ndarray:
    """Sorted-descending counts for one archetype."""
    if kind not in ARCHETYPES:
        raise ValueError(f"unknown archetype '{kind}'; expected one of {sorted(ARCHETYPES)}")
    rng = np.random.default_rng(SEED)
    parts = [lognormal_counts(rng, n, median, sigma) for n, median, sigma in ARCHETYPES[kind].values()]
    return np.sort(np.concatenate(parts))[::-1]


def n_above(counts: np.ndarray, floor: int = MIN_UMIS) -> int:
    """Barcodes the estimators would actually build a curve from."""
    return int(np.count_nonzero(counts >= floor))


def within_tolerance(rank: Optional[int], kind: str, tolerance: float = TOLERANCE_DECADES) -> Tuple[bool, float]:
    """Whether *rank* sits within *tolerance* decades of the planted boundary."""
    boundary = PLANTED_BOUNDARY[kind]
    if rank is None or not boundary:
        return False, float("nan")
    distance = abs(float(np.log10(float(rank) / float(boundary))))
    return distance <= tolerance, distance


def write(path: str, counts: np.ndarray) -> None:
    """One integer per line, as STARsolo writes UMIperCellSorted.txt."""
    with open(path, "w") as handle:
        handle.write("\n".join(str(int(c)) for c in counts) + "\n")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Write synthetic barcode-rank curves for the cell-calling checks.")
    parser.add_argument("-o", "--outdir", required=True, help="Directory the fixtures are written to.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    for kind in ARCHETYPES:
        counts = build(kind)
        path = os.path.join(args.outdir, f"{kind}_UMIperCellSorted.txt")
        write(path, counts)
        print(f"{kind}: {len(counts)} barcodes, {n_above(counts)} above {MIN_UMIS} UMIs -> {path}")

    tiny = np.asarray(TINY_COUNTS, dtype=np.int64)
    tiny_path = os.path.join(args.outdir, "tiny_UMIperCellSorted.txt")
    write(tiny_path, tiny)
    print(f"tiny: {len(tiny)} barcodes -> {tiny_path}")


if __name__ == "__main__":
    main()
