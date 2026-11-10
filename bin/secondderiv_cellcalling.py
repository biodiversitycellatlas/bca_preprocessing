#!/usr/bin/env python3
"""
Cell calling at the collapse of relative barcode weight, choosing among the bends of the
barcode-rank curve by what the barcodes on either side of them look like.

Input: the barcodes' UMI totals (``-i``, a UMIperCellSorted-style file) and, with
``--raw-dir``, the raw count matrix they were counted from.  Output: one UMI cutoff
(``-c``); barcodes whose total reaches it are cells.  The JSON (``-o``) holds the curve,
the criterion and every diagnostic of the choice.

The criterion
-------------
With ``u(r)`` the r-th largest total and ``S_r`` the cumulative sum of the top r, the
slope of log10(S_r) against log10(r) has the closed form

    g(r) = r * u(r) / S_r,

the r-th barcode's count over the mean count of the top r.  It stays near 1 on a plateau
of cells and collapses where the curve drops off a cliff.  The criterion is its rate of
change on a log-rank axis, ``dg / d log10(rank)``, differenced between adjacent windows of
constant rank ratio (WINDOW_RATIO) from prefix sums.  Bends of the curve are minima of the
criterion.  Minima are only taken from the barcodes with at least ``--min-umis`` UMIs (the
floored curve), except in tier 3 below.

Candidates and their significance
---------------------------------
Every local minimum is a candidate, with its topographic prominence.  Its significance is
measured three ways: against the criterion's robust whole-curve scale (1.4826 x MAD);
against its local noise, the median standard deviation of the criterion within
+-0.25 decades of rank over Poisson(1) bootstrap replicates of the barcodes; and by how
far the criterion climbs back on its low-count side (recovery), in local noise units.
The local scale matters because the brightest few hundred barcodes are sparse order
statistics: their wobbles look deep against the whole curve but not against their own
noise.  A minimum qualifies when it is prominent against the whole curve (>= 10) or
against its local noise with a clear recovery (>= 6 and >= 3), and in either case at
least 3 local noise scales deep.

The ambient test
----------------
The ambient profile is pooled from the barcodes with 1 to ``--ambient-max-umis`` UMIs,
and a Dirichlet-multinomial concentration is fitted to them by maximum likelihood.  As in
EmptyDrops (Lun et al., Genome Biology 2019), a barcode's statistic is the log-likelihood
of its counts under that model, against a null of AMBIENT_NULL_DRAWS simulated ambient
droplets grown one UMI at a time to AMBIENT_MAX_DEPTH.  The p-value is two-sided: a cell
whose counts sit on genes that ambient RNA is rich in is *more* probable than a typical,
overdispersed ambient draw.  Fractional (EM) counts are rounded stochastically, and
barcodes deeper than AMBIENT_MAX_DEPTH are downsampled to it.  A barcode is non-ambient
when its p-value is at most AMBIENT_P.

At a bend at rank r, ``f_above`` and ``f_below`` are the non-ambient shares of the ranks
[r / 10^0.15, r) and [r, r * 10^0.15).  A cell/empty boundary has non-ambient barcodes
above it and ambient ones below; a bend among cells has cells on both sides, and a shoulder
of debris has ambient-like barcodes on both sides, so neither shows a contrast.  A bend is
a cell boundary when ``f_above >= MIN_F_ABOVE`` and ``f_above - f_below >= MIN_CONTRAST``.

Choosing the bend
-----------------
* Tier 1: the qualifying bends that are cell boundaries.  The highest contrast wins; ties
  within CONTRAST_TIE go to the one nearest ``expected_cells`` in log rank, or without
  one to the outermost.
* Tier 2: when none is, the bends that are locally significant but fail the whole-curve
  or recovery test.  A split-pool library has no plateau of empties below its cells, so
  its boundary never climbs back; the ambient contrast stands in for the recovery test,
  and must then also have a two-proportion z of at least MIN_CONTRAST_Z.  The most
  prominent such bend wins.
* Scan: for cells that thin out with no bend at their boundary, the non-ambient share is
  read in narrow bands down the floored curve, and the boundary placed where it falls
  halfway from the cells' level to the empties'.  It is trusted only when the empties'
  level is at most SCAN_MAX_EMPTY_LEVEL, i.e. when the scan actually reaches empty
  barcodes rather than cells too shallow for the test.
* Tier 3: the minima below the floor, down to TIER3_MIN_UMIS UMIs, on tier 2's terms.
  Searching below the floor is safe only with the ambient test, which rejects a shoulder
  of debris.  It runs when nothing above placed the boundary, and also after a tier-2
  choice: when the cell types straddle the floor, a bend among them can pass tier 2
  because the test loses power on shallower cells, and a tier-3 boundary past that bend's
  own drop then wins -- the cell/empty boundary is the lowest one on the curve.  A tier-3
  minimum within the drop is its lower wall and is ignored.
* Curve rule: when the ambient test places nothing, or without ``--raw-dir``, the bend is
  chosen from the curve alone.  Of the qualifying minima within 0.1x-10x
  ``expected_cells`` the one nearest it wins, else the nearest within half a decade of
  that window, else the most prominent minimum in it.  With a matrix, minima above rank
  FALLBACK_MIN_RANK are left out of this rule.

``expected_cells`` therefore only breaks ties between equally good boundaries, and bounds
the curve rule; a warning says when the boundary found is more than 3x away from it.

Where on the bend the cutoff sits
---------------------------------
By default (``--cutoff-position inclusive``) the cutoff sits where the criterion has
climbed halfway back from the chosen minimum to its baseline, the median of the criterion
at 1.26x-2x the minimum's rank, read on the full curve so it can lie below the floor.  This
keeps the lower half of the drop, mostly cells with less RNA.  ``centroid`` cuts at the
depth-weighted centre of the minimum's basin instead.  A scanned boundary is used as is.

With ``--manual-cutoff`` the search is skipped and the given UMI threshold is written out;
the curve and its derivative are still exported.  Every random draw uses a fixed PCG64
seed, so the result depends on the inputs and nothing else.
"""

import argparse
import json
import math
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.signal import find_peaks


WINDOW_RATIO = 1.10

# Whole-curve test: prominence over the criterion's robust (MAD) noise scale.
MIN_PROMINENCE = 10.0
# Local test: prominence over the median bootstrap SD within +-0.25 decades.
MIN_LOCAL_SIGNIFICANCE = 6.0
# Floor every qualifier must clear locally, whichever test admitted it.
MIN_LOCAL_FLOOR = 3.0
# How far the criterion must climb back on the low-count side, in local noise.
MIN_RECOVERY = 3.0

LOCAL_NOISE_HALF_WIDTH = 0.25
BOOTSTRAP_REPS = 24
BOOTSTRAP_SEED = 20260904
WINDOW_EXTENSION_DECADES = 0.5

INCLUSIVE_FRACTION = 0.5
# Where the baseline the inclusive cutoff climbs back to is measured, in decades
# of rank past the minimum: 1.26x-2x its rank.
INCLUSIVE_BASELINE = (0.1, 0.3)

BASIN_FRACTION = 0.9
BASIN_SIGMAS = 1.0
MAX_RESOLUTION_DECADES = 0.05
EXPORT_POINTS = 2000
SEARCH_LO_FACTOR = 0.1
SEARCH_HI_FACTOR = 10.0

# The ambient test.  Barcodes at or below AMBIENT_MAX_UMIS make up the ambient pool;
# alpha is fitted on up to ALPHA_FIT of them.  The null is AMBIENT_NULL_DRAWS
# simulated droplets, each recorded at every depth up to AMBIENT_MAX_DEPTH.
AMBIENT_MAX_UMIS = 100
AMBIENT_PSEUDOCOUNT = 0.5
MIN_AMBIENT_POOL = 200
ALPHA_FIT = 10000
AMBIENT_NULL_DRAWS = 500
AMBIENT_MAX_DEPTH = 3000
AMBIENT_P = 0.01
# Bands either side of a bend, and how many barcodes of each are tested
BAND_DECADES = 0.15
MIN_BAND = 30
MAX_BAND = 2000
# At least this share of the barcodes above a cell boundary must be non-ambient.  Bands
# of empty droplets test at 1-6%; shallow cells just above a cliff can test as low as
# 40%, so the bar sits well between the two.
MIN_F_ABOVE = 0.2
MIN_CONTRAST = 0.05
CONTRAST_TIE = 0.05
# A second-tier bend (locally significant, but failing the whole-curve or recovery
# test) must also have a contrast this many standard errors from zero.
MIN_CONTRAST_Z = 3.0
# At most this many second-tier bends are tested, the most prominent first
TIER2_MAX_TESTS = 10
# Tier 3 looks for the boundary below the floor, down to barcodes of this many UMIs
TIER3_MIN_UMIS = 20
# When the ambient test places no boundary, the curve rule picks among the minima at or
# past this rank (the ambient bands need about as many barcodes, so none above it is tested)
FALLBACK_MIN_RANK = 100
# The scan is trusted only when its empties' level is at most this; bands of empty
# droplets test at 1-6%.
SCAN_MAX_EMPTY_LEVEL = 0.10
# A boundary this far from the expected count gets a warning that the count is off
EXPECTED_WARN_FACTOR = 3.0
AMBIENT_SEED = 20261102
# The scan used when no bend is a boundary: the non-ambient share in bands of
# +-SCAN_HALF_WIDTH decades every SCAN_STEP decades down the floored curve, each band
# subsampled to SCAN_BAND barcodes.  It needs the cells' level to stand SCAN_MIN_RANGE
# above the empties' to say anything.
SCAN_STEP_DECADES = 0.05
SCAN_HALF_WIDTH = 0.05
SCAN_BAND = 400
SCAN_MIN_RANGE = 0.2

# Guards the integer window edges against floating-point error: 11 / 1.1 is
# 10.000000000000002, and an unprotected ceil() would drop rank 11 from the curve.
_EDGE_EPS = 1e-9

# Poisson(1) CDF, long enough that 1 - CDF is below double precision.
_POISSON1_CDF = np.cumsum([math.exp(-1.0) / math.factorial(k) for k in range(25)])


def relative_weight(umis: np.ndarray) -> np.ndarray:
    """
    ``g(r) = r * u(r) / S_r`` -- the r-th barcode's count over the mean count of
    the top r.  The first derivative of log10(cumulative UMIs) with respect to
    log10(rank), in closed form, bounded in (0, 1].
    """
    ranks = np.arange(1, len(umis) + 1, dtype=np.float64)
    cum_sum = np.cumsum(umis, dtype=np.float64)
    return (ranks * umis) / cum_sum


def log_window_derivative(values: np.ndarray, ratio: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Differentiate ``values`` with respect to log10(rank) by differencing adjacent
    constant-ratio box means, ``[ceil(r/ratio), r-1]`` and ``[r, floor(r*ratio)]``,
    from prefix sums.  The denominator is the actual difference between the two
    windows' mean log10(rank).  *valid* is false where a window would be empty or
    run off the end of the curve.
    """
    n = len(values)
    ranks = np.arange(1, n + 1, dtype=np.float64)
    log_ranks = np.log10(ranks)

    csum = np.concatenate(([0.0], np.cumsum(values, dtype=np.float64)))
    lsum = np.concatenate(([0.0], np.cumsum(log_ranks, dtype=np.float64)))

    rank_i = np.arange(1, n + 1, dtype=np.int64)
    lo = np.ceil(ranks / ratio - _EDGE_EPS).astype(np.int64)
    hi = np.floor(ranks * ratio + _EDGE_EPS).astype(np.int64)
    lo = np.maximum(lo, 1)

    valid = (lo < rank_i) & (hi <= n)

    # Counts are guarded rather than masked so the vectorised gather below never
    # divides by zero on the invalid entries it is about to discard anyway.
    n_lo = np.where(valid, rank_i - lo, 1)
    n_hi = np.where(valid, hi - rank_i + 1, 1)
    lo_idx = np.clip(lo, 1, n)
    hi_idx = np.clip(hi, 1, n)

    mean_lo = (csum[rank_i - 1] - csum[lo_idx - 1]) / n_lo
    mean_hi = (csum[hi_idx] - csum[rank_i - 1]) / n_hi
    xbar_lo = (lsum[rank_i - 1] - lsum[lo_idx - 1]) / n_lo
    xbar_hi = (lsum[hi_idx] - lsum[rank_i - 1]) / n_hi

    span = xbar_hi - xbar_lo
    valid = valid & (span > 0.0)

    deriv = np.zeros(n, dtype=np.float64)
    np.divide(mean_hi - mean_lo, span, out=deriv, where=valid)
    return deriv, valid


def criterion(umis: np.ndarray, ratio: float) -> Tuple[np.ndarray, np.ndarray]:
    """``dg / d log10(rank)`` and where it is defined."""
    return log_window_derivative(relative_weight(umis.astype(np.float64)), ratio)


def robust_scale(values: np.ndarray) -> float:
    """A median absolute deviation scaled to a standard deviation."""
    if values.size == 0:
        return 0.0
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def find_candidates(deriv: np.ndarray, valid: np.ndarray) -> List[Dict[str, float]]:
    """
    Every local minimum of the criterion curve, with its topographic prominence.
    Invalid edges are filled with the curve's maximum so no basin escapes through them.
    """
    if valid.sum() < 3:
        return []

    filled = np.where(valid, deriv, np.max(deriv[valid]))
    peaks, props = find_peaks(-filled, prominence=0.0)

    candidates = []
    for i, peak in enumerate(peaks):
        if not valid[peak]:
            continue
        candidates.append(
            {
                "index": int(peak),
                "rank": int(peak) + 1,
                "prominence": float(props["prominences"][i]),
            }
        )
    return candidates


def poisson1_weights(rng: np.random.Generator, size: int) -> np.ndarray:
    """Poisson(1) draws by CDF inversion of uniforms, independent of numpy's sampler."""
    return np.searchsorted(_POISSON1_CDF, rng.random(size), side="right")


def bootstrap_noise(umis: np.ndarray, ratio: float, reps: int, seed: int, max_rank: int) -> np.ndarray:
    """
    Per-rank standard deviation of the criterion over Poisson(1)-weight bootstrap
    replicates of the barcodes, for ranks 1..max_rank.

    Each replicate repeats every barcode a Poisson(1) number of times, which keeps
    the counts sorted, so the replicate's first ranks depend only on the first
    barcodes -- truncating to the prefix the selection can reach is exact, and is
    what keeps the cost bounded on curves of millions of barcodes.  Ranks defined
    in fewer than half the replicates are NaN.
    """
    m = max(3, min(len(umis), int(max_rank)))
    prefix = umis[:m].astype(np.float64)
    rng = np.random.Generator(np.random.PCG64(seed))

    s = np.zeros(m)
    s2 = np.zeros(m)
    k = np.zeros(m)
    for _ in range(reps):
        rep = np.repeat(prefix, poisson1_weights(rng, m))
        if len(rep) < 3:
            continue
        d, v = criterion(rep, ratio)
        idx = np.flatnonzero(v[:m])
        s[idx] += d[idx]
        s2[idx] += d[idx] ** 2
        k[idx] += 1

    sd = np.full(m, np.nan)
    ok = k >= max(3, reps // 2)
    var = (s2[ok] - s[ok] ** 2 / k[ok]) / (k[ok] - 1)
    sd[ok] = np.sqrt(np.maximum(var, 0.0))
    return sd


def log_subsample(size: int, count: int) -> np.ndarray:
    """Log-spaced indices into an array of *size*."""
    if size <= count:
        return np.arange(size)
    picks = np.round(np.logspace(0.0, np.log10(size), count)).astype(np.int64) - 1
    return np.unique(np.clip(picks, 0, size - 1))


def local_noise(sd: np.ndarray, indices: Sequence[int], half_width: float) -> np.ndarray:
    """
    Median of *sd* within +-*half_width* decades of rank around each index, read on a
    log-spaced subsample so its cost does not grow with the curve.  The median, not
    the pointwise value: the knee's own location jitter inflates the pointwise SD
    on its narrow walls only.  NaN where fewer than five defined points fall in the
    window.
    """
    pos = log_subsample(len(sd), 4000)
    pos = pos[np.isfinite(sd[pos])]
    out = np.full(len(indices), np.nan)
    if pos.size == 0:
        return out
    lr = np.log10(pos + 1.0)
    vals = sd[pos]
    for j, i in enumerate(indices):
        c = math.log10(i + 1.0)
        lo, hi = np.searchsorted(lr, [c - half_width, c + half_width])
        if hi - lo >= 5:
            out[j] = float(np.median(vals[lo:hi]))
    return out


def right_recovery(deriv: np.ndarray, valid: np.ndarray, index: int) -> float:
    """
    How far the criterion climbs back above a minimum on its low-count side, before
    it either falls below the minimum again or runs out of defined ranks.  Zero for
    a minimum at the end of the defined range.
    """
    floor = deriv[index]
    best = 0.0
    n = len(deriv)
    start, step = index + 1, 64
    # Growing chunks: most walks end within a few ranks, the deepest few run the
    # whole tail, and a rank-by-rank loop would make those cost O(n) in Python.
    while start < n:
        stop = min(n, start + step)
        seg = deriv[start:stop]
        stop_here = (~valid[start:stop]) | (seg < floor)
        if stop_here.any():
            k = int(np.argmax(stop_here))
            if k > 0:
                best = max(best, float(seg[:k].max() - floor))
            return best
        best = max(best, float(seg.max() - floor))
        start, step = stop, step * 4
    return best


def basin_level(deriv: np.ndarray, index: int, prominence: float, fraction: float) -> float:
    """The contour level the basin is measured at: ``1 - fraction`` of the way out."""
    return float(deriv[index] + (1.0 - fraction) * prominence)


def contour_basin(deriv: np.ndarray, valid: np.ndarray, index: int, level: float) -> Tuple[int, int]:
    """The contiguous run around *index* that stays at or below *level*."""
    lo = hi = index
    while lo - 1 >= 0 and valid[lo - 1] and deriv[lo - 1] <= level:
        lo -= 1
    while hi + 1 < len(deriv) and valid[hi + 1] and deriv[hi + 1] <= level:
        hi += 1
    return lo, hi


def basin_centroid(deriv: np.ndarray, lo: int, hi: int, level: float) -> float:
    """Depth-weighted centroid of a basin, in log10(rank)."""
    idx = np.arange(lo, hi + 1)
    log_ranks = np.log10(idx + 1.0)
    weights = np.clip(level - deriv[idx], 0.0, None)
    total = float(weights.sum())
    if total <= 0.0:
        return float(log_ranks.mean())
    return float((weights * log_ranks).sum() / total)


def inclusive_rank(
    deriv_ext: np.ndarray,
    valid_ext: np.ndarray,
    index: int,
    centroid_rank: int,
    fraction: float,
    noise: float,
    kappa_recovery: float,
) -> Tuple[Optional[int], str]:
    """
    The first rank past the minimum at *index* where the criterion has recovered
    *fraction* of the way from the minimum to its baseline, the median of the
    criterion at 1.26x-2x the minimum's rank.

    *deriv_ext* is the criterion on every barcode with >= 1 UMI, which equals the
    floored criterion wherever both are defined (the value at rank r only involves
    barcodes ranked up to r x ratio) and keeps going where the floored curve stops,
    so a knee near the floor still has a baseline to climb back to.  The answer is
    clamped to [centroid, 2x the minimum's rank]; a bend that never recovers
    returns None.
    """
    r = index + 1
    lo = int(math.ceil(r * 10.0 ** INCLUSIVE_BASELINE[0])) - 1
    hi = int(math.floor(r * 10.0 ** INCLUSIVE_BASELINE[1])) - 1
    hi = min(hi, len(deriv_ext) - 1)
    if hi <= lo:
        return None, "truncated"
    seg = deriv_ext[lo:hi + 1][valid_ext[lo:hi + 1]]
    if seg.size < 5:
        return None, "truncated"
    base = float(np.median(seg))
    depth = base - float(deriv_ext[index])
    if not noise or depth < kappa_recovery * noise:
        return None, "no_recovery"
    level = float(deriv_ext[index]) + fraction * depth
    _, edge = contour_basin(deriv_ext, valid_ext, index, level)
    return int(np.clip(edge + 1, centroid_rank, hi + 1)), "ok"


def qualifies(c: Dict[str, float], kappa: float, kappa_local: float, kappa_recovery: float) -> bool:
    """
    The two ways a minimum can qualify, and the floor both share.

    Prominent against the whole curve, or prominent against its local noise with a
    clear climb back on its low-count side; either way at least MIN_LOCAL_FLOOR
    local noise scales deep, which is what rejects the order-statistic wobbles
    among the brightest barcodes.
    """
    if c["local_sigmas"] is None:
        return False
    whole = c["global_sigmas"] >= kappa
    local = c["local_sigmas"] >= kappa_local and c["recovery_sigmas"] >= kappa_recovery
    return (whole or local) and c["local_sigmas"] >= MIN_LOCAL_FLOOR


def _nearest(pool: Sequence[Dict[str, float]], target: float) -> Dict[str, float]:
    t = math.log10(target)
    return min(pool, key=lambda c: abs(math.log10(c["rank"]) - t))


def select_candidate(
    candidates: Sequence[Dict[str, float]],
    expected: Optional[int],
    lo_factor: float,
    hi_factor: float,
    extension: float,
) -> Tuple[Optional[Dict[str, float]], str]:
    """
    Pick the minimum that marks the cell/ambient boundary; returns ``(candidate, how)``.

    With an expected count: the qualifying minimum inside 0.1x-10x expected that
    is nearest it in log rank ("nearest_expected"); failing that, the nearest
    qualifier within *extension* decades of the window ("outside_window"); failing
    that, the most prominent minimum in the window ("shallow"), or anywhere
    ("shallow_anywhere").  Without one: the most prominent qualifier
    ("most_prominent"), else the most prominent minimum ("shallow").
    """
    if not candidates:
        return None, "no_candidates"

    qualified = [c for c in candidates if c["qualifies"]]
    if expected and expected > 0:
        lo_rank, hi_rank = expected * lo_factor, expected * hi_factor
        inner = [c for c in qualified if lo_rank <= c["rank"] <= hi_rank]
        if inner:
            return _nearest(inner, expected), "nearest_expected"
        widen = 10.0 ** extension
        outer = [c for c in qualified if lo_rank / widen <= c["rank"] <= hi_rank * widen]
        if outer:
            return _nearest(outer, expected), "outside_window"
        pool = [c for c in candidates if lo_rank <= c["rank"] <= hi_rank]
        if pool:
            return max(pool, key=lambda c: c["prominence"]), "shallow"
        return max(candidates, key=lambda c: c["prominence"]), "shallow_anywhere"

    if qualified:
        return max(qualified, key=lambda c: c["prominence"]), "most_prominent"
    return max(candidates, key=lambda c: c["prominence"]), "shallow"


ALEVIN_MATRIX = "quants_mat.mtx"


def load_ranked_matrix(raw_dir: str, matrix: Optional[str], counts: str = "SUA") -> Dict[str, object]:
    """
    A raw count directory as cells x genes CSR, rows in rank order (descending
    totals), barcodes with no UMI dropped: row ``r - 1`` is the barcode at rank r.

    Either a STARsolo raw directory (barcodes.tsv, features.tsv and *matrix*) or an
    alevin-fry quant directory (quants_mat.mtx and its row and column files), whose
    USA blocks are summed per gene over *counts* -- the same blocks the curve was
    summed over, so the two agree.  The readers are imported here, so that a run
    without ``--raw-dir`` needs nothing beyond numpy and scipy.
    """
    import scipy.sparse as sp

    if os.path.exists(os.path.join(raw_dir, ALEVIN_MATRIX)):
        import secondderiv_alevin

        mat, barcodes, columns = secondderiv_alevin.load_matrix(raw_dir)
        X, _genes = secondderiv_alevin.gene_level(mat, columns, counts)
        X = sp.csr_matrix(X, dtype=np.float64)
    else:
        import mtx_io

        name = matrix or ("UniqueAndMult-EM.mtx" if os.path.exists(os.path.join(raw_dir, "UniqueAndMult-EM.mtx"))
                          else "matrix.mtx")
        barcodes = mtx_io.read_barcodes(os.path.join(raw_dir, "barcodes.tsv"))
        n_genes = len(mtx_io.read_features(os.path.join(raw_dir, "features.tsv")))
        genes, cells, data = mtx_io._read_coordinates(os.path.join(raw_dir, name), n_genes, len(barcodes))
        X = sp.csr_matrix((data.astype(np.float64), (cells, genes)), shape=(len(barcodes), n_genes))
    totals = np.asarray(X.sum(axis=1)).ravel()
    # The curve is floored totals sorted descending; ties in it are ordered by the exact total
    order = np.lexsort((-totals, -np.floor(totals)))
    order = order[totals[order] > 0]
    return {"X": X[order], "totals": totals[order], "n_barcodes": len(barcodes)}


def integer_rows(X, rows: np.ndarray, rng: np.random.Generator, max_depth: Optional[int] = None):
    """
    Rows of *X* as whole counts: EM counts are fractional, so they are rounded
    stochastically (unbiased), and a row deeper than *max_depth* is downsampled to it
    without replacement.
    """
    sub = X[rows].tocsr(copy=True)
    sub.data = np.floor(sub.data + rng.random(sub.data.size))
    sub.eliminate_zeros()
    if max_depth:
        depth = np.asarray(sub.sum(axis=1)).ravel()
        for i in np.flatnonzero(depth > max_depth):
            s, e = sub.indptr[i], sub.indptr[i + 1]
            sub.data[s:e] = rng.multivariate_hypergeometric(sub.data[s:e].astype(np.int64), max_depth)
        sub.eliminate_zeros()
    return sub


def dm_logprob(sub, alpha: float, ap: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Each row's log-probability under DM(alpha p), and its depth; ``ap`` is alpha * p."""
    from scipy.special import gammaln

    depth = np.asarray(sub.sum(axis=1)).ravel()
    row = np.repeat(np.arange(sub.shape[0]), np.diff(sub.indptr))
    x, g = sub.data, sub.indices
    term = gammaln(x + ap[g]) - gammaln(ap[g]) - gammaln(x + 1.0)
    per_row = np.bincount(row, weights=term, minlength=sub.shape[0])
    return gammaln(depth + 1.0) + gammaln(alpha) - gammaln(depth + alpha) + per_row, depth


def fit_alpha(sub, p: np.ndarray) -> float:
    """The Dirichlet-multinomial concentration that makes the pool's barcodes most likely."""
    from scipy.optimize import minimize_scalar
    from scipy.special import gammaln

    depth = np.asarray(sub.sum(axis=1)).ravel()
    x, g = sub.data, sub.indices

    def nll(log_a):
        a = math.exp(log_a)
        ap = a * p[g]
        return -(depth.size * gammaln(a) - gammaln(depth + a).sum() + (gammaln(x + ap) - gammaln(ap)).sum())

    return float(math.exp(minimize_scalar(nll, bounds=(math.log(1e-1), math.log(1e7)), method="bounded").x))


def simulate_null(p: np.ndarray, alpha: float, draws: int, depth: int, rng: np.random.Generator) -> np.ndarray:
    """
    ``draws`` ambient droplets grown one UMI at a time: row i, column n is the DM
    log-probability of draw i after n UMIs.  Adding a UMI to a gene already holding c
    of n changes it by ``log(n+1) - log(c+1) + log(c + alpha p_g) - log(n + alpha)``,
    so each draw is one categorical sample and a cumulative sum.
    """
    ap = alpha * p
    t = np.arange(depth)
    base = np.log(t + 1.0) - np.log(t + alpha)
    null = np.zeros((draws, depth + 1))
    for i in range(draws):
        w = rng.standard_gamma(ap)
        total = w.sum()
        w = p if total <= 0 else w / total
        genes = rng.choice(len(p), size=depth, p=w)
        order = np.argsort(genes, kind="stable")
        ordered = genes[order]
        starts = np.r_[True, ordered[1:] != ordered[:-1]]
        first = np.maximum.accumulate(np.where(starts, t, 0))
        seen = np.empty(depth, dtype=np.int64)
        seen[order] = t - first
        null[i, 1:] = np.cumsum(base - np.log(seen + 1.0) + np.log(seen + ap[genes]))
    return null


def ambient_model(X, totals: np.ndarray, max_umis: float, seed: int) -> Optional[Dict[str, object]]:
    """
    The ambient profile, its Dirichlet-multinomial concentration, and the simulated
    null of the log-probability, from the barcodes with 1 to *max_umis* UMIs.  None
    when the pool is too small to estimate them.
    """
    pool = np.flatnonzero((totals >= 1) & (totals <= max_umis))
    if pool.size < MIN_AMBIENT_POOL:
        return None
    rng = np.random.Generator(np.random.PCG64(seed))
    prof = np.asarray(X[pool].sum(axis=0)).ravel() + AMBIENT_PSEUDOCOUNT
    p = prof / prof.sum()
    fit = np.sort(rng.choice(pool, size=min(ALPHA_FIT, pool.size), replace=False))
    alpha = fit_alpha(integer_rows(X, fit, rng), p)
    null = simulate_null(p, alpha, AMBIENT_NULL_DRAWS, AMBIENT_MAX_DEPTH, rng)
    return {"p": p, "alpha": alpha, "ap": alpha * p, "null": null,
            "summary": {"alpha": alpha, "pool_n": int(pool.size), "fit_n": int(fit.size), "max_umis": float(max_umis),
                        "null_draws": AMBIENT_NULL_DRAWS, "max_depth": AMBIENT_MAX_DEPTH}}


def ambient_pvalues(model: Dict[str, object], X, rows: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """
    Two-sided Monte-Carlo p-value of each barcode's counts under the ambient null, at
    its own depth.

    Two-sided: a cell whose counts sit on genes that ambient RNA is also rich in is
    *more* probable than a typical ambient draw, which with the overdispersion real
    empties show spreads its counts thinly.  A lower-tail test alone misses such cells.
    """
    sub = integer_rows(X, rows, rng, AMBIENT_MAX_DEPTH)
    lp, depth = dm_logprob(sub, model["alpha"], model["ap"])
    depth = depth.astype(np.int64)
    null = model["null"][:, depth]
    draws = model["null"].shape[0]
    lower = ((null <= lp[None, :]).sum(axis=0) + 1.0) / (draws + 1.0)
    upper = ((null >= lp[None, :]).sum(axis=0) + 1.0) / (draws + 1.0)
    pv = np.minimum(1.0, 2.0 * np.minimum(lower, upper))
    pv[depth == 0] = 1.0
    return pv


def band_contrast(model: Dict[str, object], matrix: Dict[str, object], rank: int, seed: int) -> Optional[Dict[str, float]]:
    """
    Non-ambient fractions just above and just below a bend at *rank*, each band
    subsampled to MAX_BAND barcodes.  None when either band is under MIN_BAND.
    """
    X, totals = matrix["X"], matrix["totals"]
    width = 10.0 ** BAND_DECADES
    lo = max(1, int(math.ceil(rank / width)))
    hi = min(len(totals), int(math.floor(rank * width)))
    above = np.arange(lo - 1, rank - 1)
    below = np.arange(rank - 1, hi)
    if above.size < MIN_BAND or below.size < MIN_BAND:
        return None
    rng = np.random.Generator(np.random.PCG64(seed + rank))
    out = {}
    for name, rows in (("above", above), ("below", below)):
        if rows.size > MAX_BAND:
            rows = np.sort(rng.choice(rows, size=MAX_BAND, replace=False))
        out[f"f_{name}"] = float(np.mean(ambient_pvalues(model, X, rows, rng) <= AMBIENT_P))
        out[f"n_{name}"] = int(rows.size)
    out["contrast"] = out["f_above"] - out["f_below"]
    out["z"] = contrast_z(out["f_above"], out["n_above"], out["f_below"], out["n_below"])
    return out


def contrast_z(f_above: float, n_above: int, f_below: float, n_below: int) -> float:
    """Two-proportion z of ``f_above - f_below``, with the pooled share's standard error."""
    pooled = (f_above * n_above + f_below * n_below) / (n_above + n_below)
    se = math.sqrt(pooled * (1.0 - pooled) * (1.0 / n_above + 1.0 / n_below))
    # The pooled share is 0 or 1 only when both bands are all one kind: no contrast
    return 0.0 if se == 0.0 else (f_above - f_below) / se


def select_by_ambient(
    candidates: Sequence[Dict[str, float]],
    expected: Optional[int],
    model: Dict[str, object],
    matrix: Dict[str, object],
    seed: int,
    kappa_local: float = MIN_LOCAL_SIGNIFICANCE,
) -> Tuple[Optional[Dict[str, float]], str]:
    """
    The bend that best separates non-ambient barcodes above it from ambient ones
    below it, anywhere on the floored curve, as ``(candidate, how)``.

    ``ambient_contrast`` (tier 1): a qualifying bend that is a cell boundary
    (``f_above >= MIN_F_ABOVE`` and contrast ``>= MIN_CONTRAST``); the highest contrast
    wins, ties within CONTRAST_TIE going to the one nearest the expected count, or
    without one to the outermost.

    ``ambient_tier2``: no qualifying bend is a boundary, so the bends that are locally
    significant (``local_sigmas >= kappa_local``) but failed the whole-curve or
    recovery test are tried, at most TIER2_MAX_TESTS of them by prominence.  One is a
    boundary on the same terms plus a contrast z of at least MIN_CONTRAST_Z, standing
    in for the recovery test; the most prominent wins.  Split-pool libraries need this:
    with no empty plateau below their cliff, the criterion never climbs back.

    ``(None, "no_ambient_boundary")`` when neither tier has a boundary, for the caller
    to try the scan.  Every tested bend is annotated with its band fractions and tier,
    so the JSON shows why the others lost.
    """
    def is_boundary(test: Optional[Dict[str, float]]) -> bool:
        return bool(test) and test["f_above"] >= MIN_F_ABOVE and test["contrast"] >= MIN_CONTRAST

    boundaries = []
    for c in (c for c in candidates if c["qualifies"]):
        c["ambient"] = band_contrast(model, matrix, c["rank"], seed)
        c["tier"] = 1
        if is_boundary(c["ambient"]):
            boundaries.append(c)
    if boundaries:
        best = max(c["ambient"]["contrast"] for c in boundaries)
        tied = [c for c in boundaries if c["ambient"]["contrast"] >= best - CONTRAST_TIE]
        if expected and expected > 0:
            return _nearest(tied, expected), "ambient_contrast"
        return max(tied, key=lambda c: c["rank"]), "ambient_contrast"

    pick = select_tiered([c for c in candidates if not c["qualifies"]], model, matrix, seed, kappa_local, tier=2)
    if pick is not None:
        return pick, "ambient_tier2"
    return None, "no_ambient_boundary"


def select_tiered(pool: Sequence[Dict[str, float]], model: Dict[str, object], matrix: Dict[str, object],
                  seed: int, kappa_local: float, tier: int) -> Optional[Dict[str, float]]:
    """
    The most prominent of the locally significant bends in *pool* that is a cell boundary
    with a contrast z of at least MIN_CONTRAST_Z (see ``tiered_boundaries``).
    """
    found = tiered_boundaries(pool, model, matrix, seed, kappa_local, tier)
    return found[0] if found else None


def tiered_boundaries(pool: Sequence[Dict[str, float]], model: Dict[str, object], matrix: Dict[str, object],
                      seed: int, kappa_local: float, tier: int) -> List[Dict[str, float]]:
    """
    Every locally significant bend in *pool* that is a cell boundary with a contrast z of
    at least MIN_CONTRAST_Z, most prominent first, testing at most TIER2_MAX_TESTS of
    them; each tested bend is annotated with its band fractions and *tier*.
    """
    tested = [c for c in pool if c["local_sigmas"] is not None and c["local_sigmas"] >= kappa_local]
    found = []
    for c in sorted(tested, key=lambda c: -c["prominence"])[:TIER2_MAX_TESTS]:
        c["ambient"] = band_contrast(model, matrix, c["rank"], seed)
        c["tier"] = tier
        t = c["ambient"]
        if t and t["f_above"] >= MIN_F_ABOVE and t["contrast"] >= MIN_CONTRAST and t["z"] >= MIN_CONTRAST_Z:
            found.append(c)
    return found


def bend_reach(c: Dict[str, float], deriv: np.ndarray, valid: np.ndarray, deriv_ext: np.ndarray,
               valid_ext: np.ndarray, scale: float, fraction: float, inclusive_fraction: float,
               kappa_recovery: float) -> int:
    """The last rank a floored bend's cliff covers: its inclusive cutoff, else its basin's low-count edge."""
    index = int(c["index"])
    level = basin_level(deriv, index, float(c["prominence"]), fraction)
    lo, hi = contour_basin(deriv, valid, index, level)
    centroid_rank = int(round(10.0 ** basin_centroid(deriv, lo, hi, level)))
    inc, _status = inclusive_rank(deriv_ext, valid_ext, index, centroid_rank, inclusive_fraction,
                                  c["noise"] or scale, kappa_recovery)
    return int(inc) if inc is not None else int(hi) + 1


def below_floor_candidates(curve: np.ndarray, deriv_ext: np.ndarray, valid_ext: np.ndarray, n: int, ratio: float,
                           scale: float, bootstrap_reps: int, bootstrap_seed: int) -> List[Dict[str, float]]:
    """
    The criterion's minima between the --min-umis floor and TIER3_MIN_UMIS, with their
    local noise, for tier 3.  None of them qualifies: the qualification tests stop at
    the floor, because below it a shoulder of debris can bend as sharply as a cell
    boundary.  The ambient test is what admits one -- a shoulder of debris drawn from
    ambient RNA has ambient-like barcodes on both sides.
    """
    deepest = int(np.count_nonzero(curve >= TIER3_MIN_UMIS))
    if deepest <= n:
        return []
    found = [c for c in find_candidates(deriv_ext, valid_ext) if n <= c["index"] < deepest]
    if not found:
        return []
    max_rank = min(len(curve), int(math.ceil(deepest * 10.0 ** LOCAL_NOISE_HALF_WIDTH * ratio)) + 2)
    sd = bootstrap_noise(curve, ratio, bootstrap_reps, bootstrap_seed, max_rank) if bootstrap_reps > 0 else None
    idx = [c["index"] for c in found if sd is not None and c["index"] < len(sd)]
    noise_at = dict(zip(idx, local_noise(sd, idx, LOCAL_NOISE_HALF_WIDTH))) if idx else {}
    for c in found:
        nz = noise_at.get(c["index"], math.nan)
        c["noise"] = float(nz) if np.isfinite(nz) and nz > 0 else None
        c["global_sigmas"] = c["prominence"] / scale if scale > 0 else 0.0
        c["local_sigmas"] = None if c["noise"] is None else c["prominence"] / c["noise"]
        climb = right_recovery(deriv_ext, valid_ext, c["index"])
        denom = c["noise"] or scale
        c["recovery_sigmas"] = climb / denom if denom else 0.0
        c["qualifies"] = False
        c["below_floor"] = True
    return found


def ambient_scan(model: Dict[str, object], matrix: Dict[str, object], n_floor: int, seed: int) -> Optional[Dict[str, object]]:
    """
    The boundary read off the ambient test alone, for libraries whose cells thin out
    without a bend at the boundary, where the only qualifying bends lie among the cells.

    The non-ambient share is measured in narrow bands down the floored curve and
    smoothed with a three-point running median.  The cells' level is its maximum, the
    empties' level its minimum after that, and the boundary is where it last falls
    through halfway between the two, interpolated in log rank.  ``rank`` is None when
    the two levels are within SCAN_MIN_RANGE (the curve does not show cells giving
    way to empties at all) or, with ``gate`` "empties_not_reached", when the empties'
    level is above SCAN_MAX_EMPTY_LEVEL: then the scan's lowest bands are cells too
    shallow for the test, not empty droplets, and a midpoint would fall among cells.
    """
    X, totals = matrix["X"], matrix["totals"]
    top = min(n_floor, len(totals))
    if top < 20:
        return None
    rng = np.random.Generator(np.random.PCG64(seed))
    width = 10.0 ** SCAN_HALF_WIDTH
    grid = np.unique(np.round(10.0 ** np.arange(1.0, math.log10(top) + 1e-9, SCAN_STEP_DECADES)).astype(np.int64))
    profile = []
    for r in grid:
        lo, hi = max(1, int(math.ceil(r / width))), min(len(totals), int(math.floor(r * width)))
        rows = np.arange(lo - 1, hi)
        if rows.size < MIN_BAND:
            continue
        if rows.size > SCAN_BAND:
            rows = np.sort(rng.choice(rows, size=SCAN_BAND, replace=False))
        profile.append((int(r), float(np.mean(ambient_pvalues(model, X, rows, rng) <= AMBIENT_P))))
    if len(profile) < 5:
        return None
    ranks = np.array([p[0] for p in profile], dtype=np.float64)
    f = np.array([p[1] for p in profile])
    smooth = f.copy()
    smooth[1:-1] = np.median(np.vstack([f[:-2], f[1:-1], f[2:]]), axis=0)
    peak = int(np.argmax(smooth))
    cells, empties = float(smooth[peak]), float(smooth[peak:].min())
    out = {"rank": None, "cell_level": cells, "empty_level": empties, "level": None,
           "max_empty_level": SCAN_MAX_EMPTY_LEVEL, "gate": "passed",
           "profile": [[int(r), round(float(x), 3)] for r, x in zip(ranks, f)]}
    if cells - empties < SCAN_MIN_RANGE:
        out["gate"] = "no_range"
        return out
    if empties > SCAN_MAX_EMPTY_LEVEL:
        out["gate"] = "empties_not_reached"
        return out
    level = (cells + empties) / 2.0
    k = int(np.flatnonzero(smooth[peak:] >= level).max()) + peak
    log_rank = math.log10(ranks[k])
    if k + 1 < len(smooth) and smooth[k] > smooth[k + 1]:
        t = (smooth[k] - level) / (smooth[k] - smooth[k + 1])
        log_rank += t * (math.log10(ranks[k + 1]) - math.log10(ranks[k]))
    out.update({"rank": int(np.clip(round(10.0 ** log_rank), 1, top)), "level": level})
    return out


def call_cells(
    umis: np.ndarray,
    expected_cells: Optional[int],
    *,
    extended: Optional[np.ndarray] = None,
    ratio: float = WINDOW_RATIO,
    kappa: float = MIN_PROMINENCE,
    kappa_local: float = MIN_LOCAL_SIGNIFICANCE,
    kappa_recovery: float = MIN_RECOVERY,
    bootstrap_reps: int = BOOTSTRAP_REPS,
    bootstrap_seed: int = BOOTSTRAP_SEED,
    fraction: float = BASIN_FRACTION,
    sigmas: float = BASIN_SIGMAS,
    inclusive_fraction: float = INCLUSIVE_FRACTION,
    lo_factor: float = SEARCH_LO_FACTOR,
    hi_factor: float = SEARCH_HI_FACTOR,
    extension: float = WINDOW_EXTENSION_DECADES,
    matrix: Optional[Dict[str, object]] = None,
    ambient_max_umis: float = AMBIENT_MAX_UMIS,
    ambient_seed: int = AMBIENT_SEED,
) -> Dict[str, object]:
    """
    The whole numeric path, from sorted UMI counts to both cutoff positions.

    *umis* is the floored curve (barcodes at or above --min-umis, sorted
    descending); minima are only ever taken from it.  *extended* is every barcode
    with >= 1 UMI, read only past the floor: for the recovery test, the inclusive
    cutoff and the bootstrap's local noise near the floor.  Kept out of ``main`` so
    that the comparison scripts can run this estimator beside the others.

    *matrix* is ``load_ranked_matrix``'s result for the same barcodes.  With it the
    bend is chosen by ``select_by_ambient`` over the whole floored curve; without
    it, or when no bend is a cell boundary, by the curve rule (``select_candidate``).

    ``cutoff_rank`` comes back as ``None`` when the curve is too short to
    differentiate at this window ratio.
    """
    n = len(umis)
    umis_f = umis.astype(np.float64)
    curve = umis_f if extended is None or len(extended) <= n else extended.astype(np.float64)

    deriv, valid = criterion(umis_f, ratio)
    raw_slope, raw_valid = log_window_derivative(np.log10(umis_f), ratio)

    effective_expected = None
    capped = False
    if expected_cells and expected_cells > 0:
        effective_expected = int(min(expected_cells, n))
        capped = effective_expected < expected_cells

    # The ordinary rank curve's steepest slope -- DropletUtils' inflection -- for comparison.
    raw_rank = None
    if raw_valid.any():
        searchable = raw_valid.copy()
        if effective_expected:
            ranks = np.arange(1, n + 1)
            window = searchable & (ranks >= effective_expected * lo_factor) & (ranks <= effective_expected * hi_factor)
            if window.any():
                searchable = window
        raw_rank = int(np.argmin(np.where(searchable, raw_slope, np.inf))) + 1

    result: Dict[str, object] = {
        "n": n, "deriv": deriv, "valid": valid, "raw_slope": raw_slope, "raw_valid": raw_valid,
        "raw_rank": raw_rank, "raw_umi": None if raw_rank is None else int(umis[raw_rank - 1]),
        "effective_expected": effective_expected, "capped": capped, "candidates": [], "chosen": None,
        "scale": 0.0, "how": "insufficient_curve", "prominence": None, "basin_lo": None, "basin_hi": None,
        "basin_decades": None, "resolution_lo": None, "resolution_hi": None, "resolution_decades": None,
        "cutoff_rank": None, "cutoff_umi": None, "cutoff_rank_inclusive": None, "cutoff_umi_inclusive": None,
        "inclusive_status": "not_computed", "n_qualifying": 0, "global_significance": None,
        "local_significance": None, "right_recovery": None, "local_noise": None,
        "ambient": None, "ambient_status": "no_matrix" if matrix is None else "not_computed",
        "ambient_scan": None,
    }
    if int(valid.sum()) < 3:
        return result

    scale = robust_scale(deriv[valid])
    candidates = find_candidates(deriv, valid)
    deriv_ext, valid_ext = criterion(curve, ratio) if len(curve) > n else (deriv, valid)

    # Local noise is only needed out to the deepest rank selection can reach: the
    # window widened by `extension`, capped at the floor, plus the half-width the
    # median is taken over.  The ambient test can choose a bend anywhere on the
    # floored curve, so with a matrix every bend needs it.
    reach = n if not effective_expected or matrix is not None \
        else int(min(n, math.ceil(effective_expected * hi_factor * 10.0 ** extension)))
    max_rank = min(len(curve), int(math.ceil(reach * 10.0 ** LOCAL_NOISE_HALF_WIDTH * ratio)) + 2)
    noise_at: Dict[int, float] = {}
    if candidates and bootstrap_reps > 0:
        sd = bootstrap_noise(curve, ratio, bootstrap_reps, bootstrap_seed, max_rank)
        idx = [c["index"] for c in candidates if c["index"] < len(sd)]
        if idx:
            noise_at = dict(zip(idx, local_noise(sd, idx, LOCAL_NOISE_HALF_WIDTH)))
    for c in candidates:
        nz = noise_at.get(c["index"], math.nan)
        c["noise"] = float(nz) if np.isfinite(nz) and nz > 0 else None
        c["global_sigmas"] = c["prominence"] / scale if scale > 0 else 0.0
        c["local_sigmas"] = None if c["noise"] is None else c["prominence"] / c["noise"]
        climb = right_recovery(deriv_ext, valid_ext, c["index"])
        denom = c["noise"] or scale
        c["recovery_sigmas"] = climb / denom if denom else 0.0
        c["qualifies"] = qualifies(c, kappa, kappa_local, kappa_recovery)

    chosen, how, scan = None, None, None
    if matrix is not None:
        model = ambient_model(matrix["X"], matrix["totals"], ambient_max_umis, ambient_seed)
        if model is None:
            result["ambient_status"] = "pool_too_small"
        else:
            result["ambient"] = model["summary"]
            chosen, how = select_by_ambient(candidates, effective_expected, model, matrix, ambient_seed,
                                            kappa_local=kappa_local)
            if chosen is not None:
                result["ambient_status"] = "ok"
            else:
                # No bend is a boundary: read the boundary off the ambient test instead
                scan = ambient_scan(model, matrix, n, ambient_seed)
                result["ambient_scan"] = scan
                if scan and scan["rank"]:
                    how = result["ambient_status"] = "ambient_scan"
                elif scan and scan["gate"] == "empties_not_reached":
                    result["ambient_status"] = "scan_rejected"
                else:
                    result["ambient_status"] = how
            if (chosen is None and how != "ambient_scan") or how == "ambient_tier2":
                # Tier 3: the boundary may lie below the floor, as in a shallow split-pool library.
                # After a tier-2 pick too: the cell/empty boundary is the lowest one, and a split-pool
                # library whose cells straddle the floor has a bend among them (deep types above,
                # shallow ones below) that passes tier 2 because the test loses power with depth.
                below = below_floor_candidates(curve, deriv_ext, valid_ext, n, ratio, scale, bootstrap_reps,
                                               bootstrap_seed)
                candidates.extend(below)
                found = tiered_boundaries(below, model, matrix, ambient_seed, kappa_local, tier=3)
                if how == "ambient_tier2":
                    # Only a boundary past the tier-2 bend's own cliff counts; one inside it is the
                    # cliff's lower wall
                    reach = bend_reach(chosen, deriv, valid, deriv_ext, valid_ext, scale, fraction,
                                       inclusive_fraction, kappa_recovery)
                    found = [c for c in found if c["rank"] > reach]
                if found:
                    chosen, how = found[0], "ambient_tier3"
                    result["ambient_status"] = "ok"
    if chosen is None and how != "ambient_scan":
        # The curve rule, on the floored curve's minima only.  After the ambient test has
        # failed, the brightest FALLBACK_MIN_RANK barcodes are left out too: their sparse
        # order statistics make prominent minima that a low or missing expected count
        # would otherwise land on.
        floored = [c for c in candidates if not c.get("below_floor")]
        if matrix is not None:
            floored = [c for c in floored if c["rank"] >= FALLBACK_MIN_RANK] or floored
        chosen, how = select_candidate(floored, effective_expected, lo_factor, hi_factor, extension)

    # A tier-3 bend lies below the floor: its basin is read on the extended curve
    below_floor = how != "ambient_scan" and chosen is not None and int(chosen["index"]) >= n
    d_use, v_use = (deriv_ext, valid_ext) if below_floor else (deriv, valid)
    if how == "ambient_scan":
        index = int(scan["rank"]) - 1
        basin_lo = basin_hi = index
        centroid_logx = float(np.log10(index + 1.0))
        basin_decades = None
        prominence = None
    elif chosen is None:
        index = int(np.argmin(np.where(valid, deriv, np.inf)))
        basin_lo = basin_hi = index
        centroid_logx = float(np.log10(index + 1.0))
        basin_decades = None
        prominence = None
        how = "global_argmin"
    else:
        index = int(chosen["index"])
        prominence = float(chosen["prominence"])
        level = basin_level(d_use, index, prominence, fraction)
        basin_lo, basin_hi = contour_basin(d_use, v_use, index, level)
        centroid_logx = basin_centroid(d_use, basin_lo, basin_hi, level)
        basin_decades = float(np.log10((basin_hi + 1.0) / (basin_lo + 1.0)))

    # Resolution: the ranks within one whole-curve noise scale of the minimum, i.e.
    # the ranks not distinguishable from the chosen one.
    # A scanned boundary is not a bend, so it has no resolution on the criterion.
    resolution_lo = resolution_hi = index
    resolution_decades = None
    if scale > 0.0 and how != "ambient_scan":
        resolution_lo, resolution_hi = contour_basin(d_use, v_use, index, float(d_use[index] + sigmas * scale))
        resolution_decades = float(np.log10((resolution_hi + 1.0) / (resolution_lo + 1.0)))

    cutoff_rank = int(np.clip(round(10.0 ** centroid_logx), 1, len(curve) if below_floor else n))

    inc_rank, inc_status = None, "no_minimum"
    if how == "ambient_scan":
        # The scan already places the cut where cells give way to empties
        inc_rank, inc_status = cutoff_rank, "ambient_scan"
    elif chosen is not None:
        inc_rank, inc_status = inclusive_rank(
            deriv_ext, valid_ext, index, cutoff_rank, inclusive_fraction,
            chosen["noise"] or scale, kappa_recovery,
        )
    inc_umi = None if inc_rank is None else int(curve[min(inc_rank, len(curve)) - 1])

    result.update(
        {
            "scale": scale, "candidates": candidates, "chosen": chosen, "how": how, "prominence": prominence,
            "basin_lo": basin_lo, "basin_hi": basin_hi, "basin_decades": basin_decades,
            "resolution_lo": resolution_lo, "resolution_hi": resolution_hi, "resolution_decades": resolution_decades,
            "cutoff_rank": cutoff_rank, "cutoff_umi": int(curve[cutoff_rank - 1]),
            "cutoff_rank_inclusive": inc_rank, "cutoff_umi_inclusive": inc_umi, "inclusive_status": inc_status,
            "n_qualifying": sum(1 for c in candidates if c["qualifies"]),
            "global_significance": None if chosen is None else chosen["global_sigmas"],
            "local_significance": None if chosen is None else chosen["local_sigmas"],
            "right_recovery": None if chosen is None else chosen["recovery_sigmas"],
            "local_noise": None if chosen is None else chosen["noise"],
        }
    )
    return result


def write_fallback(sample_id: str, out_json: str, out_cutoff: str, reason: str, cutoff: int = 0) -> None:
    """Write valid fallback outputs for datasets that are too small to analyze."""
    export_data = {
        sample_id: {
            "logX": [], "logY": [], "customdata": [], "derivX": [], "derivY": [],
            "threshold_logX": None, "threshold_logY": None, "threshold_umi": cutoff,
            "status": "warning", "message": reason,
        }
    }
    with open(out_json, "w") as f:
        json.dump(export_data, f, separators=(",", ":"))
    with open(out_cutoff, "w") as f:
        f.write(str(cutoff))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Call cells at the collapse of relative barcode weight, choosing the bend whose two sides "
                    "differ most in how ambient-like their barcodes are (with --raw-dir), and export JSON."
    )
    parser.add_argument("-i", "--input", required=True, help="Input UMIperCellSorted.txt file.")
    parser.add_argument("-s", "--sample-id", required=True, help="Sample ID used as the JSON top-level key.")
    parser.add_argument("-o", "--out-json", default="knee_data.json", help="Output JSON file.")
    parser.add_argument("-c", "--out-cutoff", default="cutoff.txt", help="Output cutoff text file.")
    parser.add_argument("-e", "--expected-cells", type=int, default=None,
                        help="Expected number of cells. With --raw-dir it only breaks ties between equally "
                             "good cell boundaries; when the bend is chosen from the curve alone, the nearest "
                             f"qualifying minimum inside {SEARCH_LO_FACTOR:g}x-{SEARCH_HI_FACTOR:g}x this value wins.")
    parser.add_argument("--raw-dir", default=None,
                        help="Raw count directory the curve was counted from: a STARsolo raw directory "
                             "(barcodes.tsv, features.tsv and the matrix) or an alevin-fry quant directory "
                             "(quants_mat.mtx and its row and column files). Enables the ambient test that "
                             "chooses the bend; without it the bend is chosen from the curve alone.")
    parser.add_argument("--matrix", default=None,
                        help="STARsolo matrix in --raw-dir. Default: UniqueAndMult-EM.mtx when present, else "
                             "matrix.mtx. It must be the matrix the -i curve comes from.")
    parser.add_argument("--counts", default="SUA", choices=["SUA", "SA", "S", "UA", "U"],
                        help="alevin-fry --raw-dir only: the USA blocks summed per gene, the same blocks the -i "
                             "curve was summed over (default: SUA).")
    parser.add_argument("--ambient-max-umis", type=float, default=AMBIENT_MAX_UMIS,
                        help="Barcodes with at most this many UMIs make up the ambient pool.")
    parser.add_argument("-m", "--manual-cutoff", type=int, default=None,
                        help="Skip the search and write this UMI cutoff instead. The curve and its "
                             "derivative are still exported, with the threshold marked at the given value.")
    parser.add_argument("--cutoff-position", choices=["inclusive", "centroid"], default="inclusive",
                        help="Where on the chosen bend the cutoff sits: 'inclusive' (default) is the low-count "
                             "end of the drop, keeping the cells with less RNA on its lower half; 'centroid' "
                             "is the middle of the drop.")
    parser.add_argument("--inclusive-fraction", type=float, default=INCLUSIVE_FRACTION,
                        help="How far the criterion must have climbed back from the minimum towards its "
                             "baseline for the inclusive cutoff (0-1).")
    parser.add_argument("--window-ratio", type=float, default=WINDOW_RATIO,
                        help="Ratio of the box-mean windows the derivative is taken over.")
    parser.add_argument("--min-prominence", type=float, default=MIN_PROMINENCE,
                        help="Whole-curve test: prominence in robust noise scales.")
    parser.add_argument("--min-local-significance", type=float, default=MIN_LOCAL_SIGNIFICANCE,
                        help="Local test: prominence in local bootstrap noise scales.")
    parser.add_argument("--min-recovery", type=float, default=MIN_RECOVERY,
                        help="Local test: climb back on the low-count side, in local noise scales.")
    parser.add_argument("--bootstrap-reps", type=int, default=BOOTSTRAP_REPS,
                        help="Bootstrap replicates for the local noise scale.")
    parser.add_argument("--bootstrap-seed", type=int, default=BOOTSTRAP_SEED, help="Seed for the bootstrap.")
    parser.add_argument("--export-points", type=int, default=EXPORT_POINTS,
                        help="Points of the derivative curve written to the JSON for plotting.")
    parser.add_argument("--min-umis", type=int, default=100, help="Minimum UMI threshold used to retain entries before analysis.")
    parser.add_argument("--min-points", type=int, default=30, help="Minimum number of retained points required for an unflagged result.")
    parser.add_argument("--fallback-cutoff", type=int, default=0, help="Cutoff written when the dataset is too small to analyze.")
    return parser.parse_args()


def describe(called: Dict[str, object], args: argparse.Namespace, umis: np.ndarray, n: int) -> List[str]:
    """Everything the caller should be told about how the cutoff was reached."""
    warnings: List[str] = []
    how = called["how"]
    chosen = called["chosen"]
    if n < args.min_points:
        warnings.append(f"Only {n} barcodes are above {args.min_umis} UMIs; interpret the cutoff cautiously.")
    if called["capped"]:
        warnings.append(
            f"Only {n} barcodes are above {args.min_umis} UMIs, far fewer than the {args.expected_cells} "
            f"expected cells, so the window was built from those {n} instead. This sample has nowhere "
            "near as many cells as expected; treat the cutoff, and the sample, with suspicion."
        )
    status = called.get("ambient_status")
    scan = called.get("ambient_scan")
    if status == "ambient_scan":
        warnings.append(
            f"No bend separates cells from empty droplets, so the cutoff is where the share of barcodes unlike "
            f"ambient RNA falls halfway from the cells' {scan['cell_level']:.0%} to the empties' "
            f"{scan['empty_level']:.0%}, at rank {scan['rank']}. The cells thin out gradually here, so "
            "any single cutoff is a compromise; check it against the barcode-rank plot."
        )
    elif status not in (None, "ok", "no_matrix"):
        best = max((c["ambient"]["contrast"] for c in called["candidates"] if c.get("ambient")), default=None)
        reason = {
            "pool_too_small": f"too few barcodes with at most {args.ambient_max_umis:g} UMIs to model ambient RNA",
            "no_ambient_boundary": "no bend separates non-ambient barcodes above it from ambient ones below it"
                                   + ("" if best is None else f" (best contrast {best:.2f} of "
                                                              f"{MIN_CONTRAST:g} needed)")
                                   + ("" if not scan else f", and along the curve the non-ambient share only "
                                                         f"ranges {scan['empty_level']:.0%}-{scan['cell_level']:.0%}"),
            "scan_rejected": "no bend separates non-ambient barcodes above it from ambient ones below it"
                             + ("" if best is None else f" (best contrast {best:.2f} of {MIN_CONTRAST:g} needed)")
                             + ("" if not scan else
                                f", and the scan's lowest bands are still {scan['empty_level']:.0%} unlike ambient "
                                f"RNA, above the {SCAN_MAX_EMPTY_LEVEL:.0%} of empty droplets: they are cells too "
                                "shallow for the test, so the scan's midpoint would fall among cells"),
        }.get(status, status)
        warnings.append(f"The ambient test could not place the boundary ({reason}), so the bend was chosen "
                        "from the curve alone.")
    if how in ("ambient_contrast", "ambient_tier2", "ambient_tier3"):
        t = chosen["ambient"]
        if how == "ambient_tier2":
            warnings.append(
                f"No qualifying bend is a cell boundary, so a locally significant one that failed the "
                f"whole-curve or recovery test was taken (rank {chosen['rank']}, {int(umis[chosen['rank'] - 1])} "
                f"UMIs; local {chosen['local_sigmas']:.1f} sigma, recovery {chosen['recovery_sigmas']:.1f}): "
                f"{t['f_above']:.0%} of the barcodes above it and {t['f_below']:.0%} below are unlike ambient "
                f"RNA (contrast {t['contrast']:.2f}, z {t['z']:.1f}). Libraries without an empty-droplet "
                "plateau, such as split-pool ones, look like this."
            )
        elif how == "ambient_tier3":
            warnings.append(
                f"No bend above {args.min_umis} UMIs is a cell boundary, so one below that floor was taken "
                f"(rank {chosen['rank']}, {int(umis[chosen['rank'] - 1])} UMIs; local "
                f"{chosen['local_sigmas']:.1f} sigma): {t['f_above']:.0%} of the barcodes above it and "
                f"{t['f_below']:.0%} below are unlike ambient RNA (contrast {t['contrast']:.2f}, z {t['z']:.1f}). "
                "The cells of this library are shallow; check the cutoff against the barcode-rank plot."
            )
        if args.expected_cells and not (1 / EXPECTED_WARN_FACTOR <= chosen["rank"] / args.expected_cells
                                        <= EXPECTED_WARN_FACTOR):
            warnings.append(
                f"The cell boundary is at rank {chosen['rank']}, more than {EXPECTED_WARN_FACTOR:g}x from the "
                f"{args.expected_cells} expected cells: {t['f_above']:.0%} of the barcodes above it and "
                f"{t['f_below']:.0%} below it are unlike ambient RNA. expected_cells is probably off."
            )
    elif how == "ambient_scan":
        if args.expected_cells and not (1 / EXPECTED_WARN_FACTOR <= scan["rank"] / args.expected_cells
                                        <= EXPECTED_WARN_FACTOR):
            warnings.append(
                f"The boundary is at rank {scan['rank']}, more than {EXPECTED_WARN_FACTOR:g}x from the "
                f"{args.expected_cells} expected cells. expected_cells is probably off."
            )
    elif not args.expected_cells:
        warnings.append(
            "No expected cell count given, so the most prominent qualifying minimum on the whole curve "
            "was taken; a curve with several bends may need one."
        )
    if how in ("shallow", "shallow_anywhere"):
        warnings.append(
            "No minimum passed the qualification tests, so the most prominent of a shallow field was "
            "taken. There may be no defensible knee in this library; check the barcode-rank plot, or set a "
            "manual cutoff."
        )
    elif how == "outside_window":
        warnings.append(
            f"No qualifying minimum lies within {SEARCH_LO_FACTOR:g}-{SEARCH_HI_FACTOR:g}x the "
            f"{args.expected_cells} expected cells; the nearest one just outside, at rank {chosen['rank']}, "
            "was taken. The expected cell count is probably off by more than tenfold."
        )
    elif how == "global_argmin":
        warnings.append("The criterion curve has no local minimum, so its lowest point was taken.")
    if chosen is not None and chosen["qualifies"] and chosen["global_sigmas"] < args.min_prominence:
        warnings.append(
            f"The chosen bend is shallow against the whole curve ({chosen['global_sigmas']:.1f} of "
            f"{args.min_prominence:g} noise scales) and was accepted on its local noise "
            f"({chosen['local_sigmas']:.1f}) and recovery ({chosen['recovery_sigmas']:.1f}), so the cutoff "
            "is only loosely determined."
        )
    if chosen is not None and how == "nearest_expected":
        qualifying = [c for c in called["candidates"] if c["qualifies"]]
        deepest = max(qualifying, key=lambda c: c["prominence"]) if qualifying else None
        if deepest is not None and deepest["index"] != chosen["index"]:
            warnings.append(
                f"The expected cell count moved the answer: rank {chosen['rank']} was chosen for being "
                f"nearest {called['effective_expected']}, over the more prominent qualifying minimum at "
                f"rank {deepest['rank']} ({int(umis[deepest['rank'] - 1])} UMIs)."
            )
    resolution = called["resolution_decades"]
    if resolution is not None and resolution > MAX_RESOLUTION_DECADES:
        warnings.append(
            f"Ranks {int(called['resolution_lo']) + 1}-{int(called['resolution_hi']) + 1} "
            f"({resolution:.2f} decades) are indistinguishable from the chosen one at this curve's noise "
            "level, so the cutoff is only loosely determined."
        )
    return warnings


def main() -> None:
    args = parse_args()

    try:
        data = np.loadtxt(args.input, dtype=int)
        data = np.atleast_1d(data).flatten()

        all_counts = np.sort(data)[::-1]
        all_counts = all_counts[all_counts > 0]
        umis = all_counts[all_counts >= args.min_umis]
        n = len(umis)

        if n < 10:
            msg = f"Too few points above {args.min_umis} UMIs for reliable cutoff detection ({n} found)."
            print(f"Warning: {msg}", file=sys.stderr)
            # A manually chosen cutoff stands even when the curve is too short to
            # analyse; there is nothing to detect, only a threshold to apply.
            fallback = args.fallback_cutoff if args.manual_cutoff is None else args.manual_cutoff
            write_fallback(args.sample_id, args.out_json, args.out_cutoff, msg, cutoff=fallback)
            return

        if args.window_ratio <= 1.0:
            raise ValueError(f"--window-ratio must be greater than 1 (got {args.window_ratio}).")
        if not 0.0 < args.inclusive_fraction < 1.0:
            raise ValueError(f"--inclusive-fraction must lie between 0 and 1 (got {args.inclusive_fraction}).")

        ranks = np.arange(1, n + 1)
        cum_sum = np.cumsum(umis.astype(np.float64))
        log_x, log_y = np.log10(ranks.astype(np.float64)), np.log10(cum_sum)

        warnings: List[str] = []
        matrix = None
        if args.raw_dir and args.manual_cutoff is None:
            if os.path.exists(os.path.join(args.raw_dir, ALEVIN_MATRIX)):
                name = f"{ALEVIN_MATRIX} ({args.counts})"
            else:
                name = args.matrix or ("UniqueAndMult-EM.mtx" if os.path.exists(
                    os.path.join(args.raw_dir, "UniqueAndMult-EM.mtx")) else "matrix.mtx")
            matrix = load_ranked_matrix(args.raw_dir, args.matrix, args.counts)
            # The bands are read by rank, so the matrix must give the curve's ranks
            floored = np.floor(matrix["totals"]).astype(np.int64)
            floored = floored[floored > 0]
            m = min(len(floored), len(all_counts))
            mismatch = float(np.mean(np.abs(floored[:m] - all_counts[:m]) > 1)) if m else 1.0
            if mismatch > 0.01:
                warnings.append(
                    f"The curve does not match the totals of {name} in {args.raw_dir} ({mismatch:.1%} of ranks "
                    "differ), so the ambient test was skipped and the bend chosen from the curve alone. Give "
                    "the matrix the curve was counted from (--matrix, or --counts for alevin-fry)."
                )
                matrix = None

        called = call_cells(
            umis, args.expected_cells, extended=all_counts, ratio=args.window_ratio,
            kappa=args.min_prominence, kappa_local=args.min_local_significance,
            kappa_recovery=args.min_recovery, bootstrap_reps=args.bootstrap_reps,
            bootstrap_seed=args.bootstrap_seed, inclusive_fraction=args.inclusive_fraction,
            matrix=matrix, ambient_max_umis=args.ambient_max_umis,
        )
        deriv, valid = called["deriv"], called["valid"]
        candidates = called["candidates"]
        how = called["how"]
        position = args.cutoff_position

        if args.manual_cutoff is not None:
            how = position = "manual"
            final_cutoff = int(args.manual_cutoff)
            kept = int(np.count_nonzero(umis >= final_cutoff))
            cutoff_rank = int(np.clip(kept, 1, n))
            if kept == 0:
                warnings.append(
                    f"Manual cutoff of {final_cutoff} UMIs is above every barcode in this sample "
                    f"(highest is {int(umis[0])}); no cells would be retained."
                )
            elif final_cutoff < args.min_umis:
                # The plotted curve stops at min_umis, so a lower threshold falls past its
                # end and the marker is clamped there; the cutoff itself is still applied
                # to every barcode in the matrix.
                warnings.append(
                    f"Manual cutoff of {final_cutoff} UMIs is below the {args.min_umis} UMIs this curve "
                    "is plotted down to, so the marker is clamped to its tail; the cutoff is still "
                    "applied to every barcode."
                )
            status = "warning" if warnings else "ok"
            message = " ".join(warnings) if warnings else (
                f"Manual cutoff of {final_cutoff} UMIs applied; the search was skipped and {kept} "
                "barcodes are retained."
            )

        elif called["cutoff_rank"] is None:
            msg = (
                f"Only {int(valid.sum())} of {n} barcodes above {args.min_umis} UMIs have a defined "
                f"derivative at --window-ratio {args.window_ratio:g}; the curve is too short to differentiate."
            )
            print(f"Warning: {msg}", file=sys.stderr)
            write_fallback(args.sample_id, args.out_json, args.out_cutoff, msg, cutoff=args.fallback_cutoff)
            return

        else:
            cutoff_rank = int(called["cutoff_rank"])
            final_cutoff = int(called["cutoff_umi"])
            if position == "inclusive":
                if called["cutoff_rank_inclusive"] is None:
                    warnings.append(
                        f"The inclusive cutoff could not be placed ({called['inclusive_status']}: the bend "
                        "does not climb back clearly on its low-count side), so the centroid was used."
                    )
                    position = "centroid"
                else:
                    final_cutoff = int(called["cutoff_umi_inclusive"])
                    cutoff_rank = int(called["cutoff_rank_inclusive"])
            if cutoff_rank > n:
                warnings.append(
                    f"The cutoff ({final_cutoff} UMIs) lies below the {args.min_umis} UMIs this curve is "
                    "plotted down to. Its marker is clamped to the tail of the plot; the cutoff is still "
                    "applied to every barcode."
                )
                cutoff_rank = n
            warnings.extend(describe(called, args, all_counts, n))
            status = "warning" if warnings else "ok"
            message = " ".join(warnings) if warnings else "Cutoff calculated successfully."

        export_idx = log_subsample(n, args.export_points)
        export_idx = export_idx[valid[export_idx]]
        deriv_ceiling = float(np.max(deriv[valid])) if valid.any() else 0.0
        basin_lo, basin_hi = called["basin_lo"], called["basin_hi"]

        def opt(v, cast=float):
            return None if v is None else cast(v)

        entry = {
            "logX": log_x.tolist(),
            "logY": log_y.tolist(),
            "customdata": [[int(r), int(u)] for r, u in zip(ranks, umis)],
            "derivX": log_x[export_idx].tolist(),
            "derivY": (deriv[export_idx] - deriv_ceiling).tolist(),
            "threshold_logX": float(log_x[cutoff_rank - 1]),
            "threshold_logY": float(log_y[cutoff_rank - 1]),
            "threshold_umi": final_cutoff,
            "status": status,
            "message": message,
            # Additive keys; consumers that do not know them ignore them.
            "criterion": "dg_dlogrank",
            "script_version": 6,
            "selection_tier": (1 if how == "ambient_contrast" else 2 if how == "ambient_tier2"
                               else 3 if how == "ambient_tier3" else "scan" if how == "ambient_scan"
                               else None if how == "manual" else "curve"),
            "ambient_status": called["ambient_status"],
            "ambient_model": called["ambient"],
            "ambient_scan": called["ambient_scan"],
            "ambient_tests": [
                {"rank": c["rank"], "umi": int(all_counts[c["rank"] - 1]), "tier": c.get("tier"), **c["ambient"]}
                for c in sorted(candidates, key=lambda c: c["rank"]) if c.get("ambient")
            ],
            "selection_rule": how,
            "cutoff_position": position,
            "cutoff_rank": cutoff_rank,
            "cutoff_rank_centroid": opt(called["cutoff_rank"], int),
            "cutoff_umi_centroid": opt(called["cutoff_umi"], int),
            "cutoff_rank_inclusive": opt(called["cutoff_rank_inclusive"], int),
            "cutoff_umi_inclusive": opt(called["cutoff_umi_inclusive"], int),
            "inclusive_fraction": float(args.inclusive_fraction),
            "inclusive_status": called["inclusive_status"],
            "window_ratio": float(args.window_ratio),
            "expected_cells_effective": called["effective_expected"],
            "bootstrap_reps": int(args.bootstrap_reps),
            "bootstrap_seed": int(args.bootstrap_seed),
            "global_significance": opt(called["global_significance"]),
            "local_significance": opt(called["local_significance"]),
            "local_noise": opt(called["local_noise"]),
            "right_recovery": opt(called["right_recovery"]),
            "n_qualifying": int(called["n_qualifying"]),
            "prominence": called["prominence"],
            "basin_rank_lo": None if basin_lo is None else int(basin_lo) + 1,
            "basin_rank_hi": None if basin_hi is None else int(basin_hi) + 1,
            "basin_width_decades": called["basin_decades"],
            "resolution_rank_lo": None if called["resolution_lo"] is None else int(called["resolution_lo"]) + 1,
            "resolution_rank_hi": None if called["resolution_hi"] is None else int(called["resolution_hi"]) + 1,
            "resolution_decades": called["resolution_decades"],
            "n_candidates": len(candidates),
            # Only the deepest few are worth shipping; the count above is the truth.
            "candidates": [
                {
                    "rank": c["rank"],
                    "prominence": round(c["prominence"], 6),
                    "umi": int(all_counts[c["rank"] - 1]),
                    "global_sigmas": round(c["global_sigmas"], 3),
                    "local_sigmas": None if c["local_sigmas"] is None else round(c["local_sigmas"], 3),
                    "recovery_sigmas": round(c["recovery_sigmas"], 3),
                    "qualifies": bool(c["qualifies"]),
                    "contrast": None if not c.get("ambient") else round(c["ambient"]["contrast"], 3),
                    "contrast_z": None if not c.get("ambient") else round(c["ambient"]["z"], 2),
                    "tier": c.get("tier"),
                }
                for c in sorted(candidates, key=lambda c: -c["prominence"])[:20]
            ],
            "raw_slope_derivY": called["raw_slope"][export_idx].tolist(),
            "raw_slope_rank": called["raw_rank"],
            "raw_slope_umi": called["raw_umi"],
        }

        with open(args.out_json, "w") as f:
            json.dump({args.sample_id: entry}, f, separators=(",", ":"))

        with open(args.out_cutoff, "w") as f:
            f.write(str(final_cutoff))

        if args.manual_cutoff is not None:
            detail = "manual cutoff"
        else:
            detail = (f"rule={how}, position={position}, {called['n_qualifying']} of {len(candidates)} "
                      "minima qualify")
            if called["effective_expected"] is not None:
                detail += f", expected_cells={args.expected_cells}"
                if called["capped"]:
                    detail += f" capped to {called['effective_expected']}"
        print(
            f"Processed {args.sample_id} -> UMI Cutoff: {final_cutoff} "
            f"(rank {cutoff_rank} of {n} barcodes above {args.min_umis} UMIs; {detail})",
            file=sys.stderr,
        )

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
