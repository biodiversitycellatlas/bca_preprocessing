#!/usr/bin/env python3
"""
Pre-calculate a second-derivative-based cell-calling cutoff.

Cells are separated from empty droplets at the sharpest downward bend of the
log10(cumulative UMIs) vs log10(barcode rank) curve.  Two things keep that bend
findable on real data:

* the curve is resampled onto a grid that is uniform in log10(rank) before it is
  differentiated, so the derivative is not divided by a step that shrinks as
  ~1/rank and amplifies tail noise in proportion to rank; and
* the minimum is searched for within a rank window around the expected cell
  count, because the cumulative curve of a sample with a deep ambient tail keeps
  bending well past the real cells and an unbounded search settles out there.
"""

import argparse
import json
import sys
from typing import Optional, Tuple

import numpy as np


GRID_POINTS = 2000
SMOOTH_FRACTION = 0.02
SEARCH_LO_FACTOR = 0.1
SEARCH_HI_FACTOR = 10.0


def moving_average(data: np.ndarray, window_size: int) -> np.ndarray:
    """
    Smooth the array using a centered moving average.

    The edges are replicated rather than zero-padded: ``np.convolve(..., mode="same")``
    averages in zeros at both ends, which drags the curve towards zero there and
    plants a spurious minimum that the inflection search then locks onto.
    """
    if window_size < 3:
        return data
    pad = window_size // 2
    padded = np.pad(data, pad, mode="edge")
    return np.convolve(padded, np.ones(window_size) / window_size, mode="valid")


def resample_log_rank(
    log_x: np.ndarray, log_y: np.ndarray, n_grid: int, smooth_window: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Resample the curve onto a grid that is uniform in ``log10(rank)`` and smooth it.

    The raw curve is indexed by rank, so the spacing of ``log10(rank)`` collapses as
    ~1/rank towards the tail.  Differentiating against that spacing divides by an
    ever smaller step, amplifying the integer quantisation of the UMI counts by a
    factor that grows linearly with rank -- which buries the knee under tail noise
    and sends the minimum of the second derivative into the ambient droplets.
    Resampling onto a uniform grid removes that amplification.
    """
    grid_x = np.linspace(log_x[0], log_x[-1], n_grid)
    grid_y = np.interp(grid_x, log_x, log_y)
    return grid_x, moving_average(grid_y, smooth_window)


def calculate_derivatives(grid_x: np.ndarray, grid_y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute first and second derivatives on the uniform ``log10(rank)`` grid.
    """
    d1 = np.gradient(grid_y, grid_x)
    d2 = np.gradient(d1, grid_x)
    return d1, d2


def detect_inflection_point(
    second_deriv: np.ndarray,
    ranks: np.ndarray,
    trim_edge: int,
    expected_cells: Optional[int],
    lo_factor: float,
    hi_factor: float,
) -> Tuple[int, bool]:
    """
    Index of the most negative second-derivative point, i.e. the sharpest downward
    bend of the cumulative curve.

    Returns ``(index, bounded)``, where *bounded* records whether the search was
    restricted to a rank window around *expected_cells*.  Without that window the
    global minimum is not reliably the knee: on a sample with a deep ambient tail
    the cumulative curve keeps bending well past the real cells, and the search
    walks out into the tail and returns a cutoff that retains every barcode.
    """
    searchable = np.ones(len(second_deriv), dtype=bool)
    if len(second_deriv) > (trim_edge * 2):
        searchable[:trim_edge] = False
        searchable[-trim_edge:] = False

    bounded = False
    if expected_cells and expected_cells > 0:
        window = (
            searchable
            & (ranks >= expected_cells * lo_factor)
            & (ranks <= expected_cells * hi_factor)
        )
        if window.any():
            searchable, bounded = window, True

    if not searchable.any():
        searchable = np.ones(len(second_deriv), dtype=bool)

    return int(np.argmin(np.where(searchable, second_deriv, np.inf))), bounded


def write_fallback(sample_id: str, out_json: str, out_cutoff: str, reason: str, cutoff: int = 0) -> None:
    """
    Write valid fallback outputs for datasets that are too small to analyze.
    """
    export_data = {
        sample_id: {
            "logX": [],
            "logY": [],
            "customdata": [],
            "derivX": [],
            "derivY": [],
            "threshold_logX": None,
            "threshold_logY": None,
            "threshold_umi": cutoff,
            "status": "warning",
            "message": reason,
        }
    }

    with open(out_json, "w") as f:
        json.dump(export_data, f, separators=(",", ":"))

    with open(out_cutoff, "w") as f:
        f.write(str(cutoff))


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Calculate a second-derivative cell-calling cutoff and export JSON.")
    parser.add_argument("-i", "--input", required=True, help="Input UMIperCellSorted.txt file.")
    parser.add_argument("-s", "--sample-id", required=True, help="Sample ID used as the JSON top-level key.")
    parser.add_argument( "-o", "--out-json", default="knee_data.json", help="Output JSON file.")
    parser.add_argument("-c", "--out-cutoff", default="cutoff.txt", help="Output cutoff text file.")
    parser.add_argument("-e", "--expected-cells", type=int, default=None,
                        help="Expected number of cells. The knee is searched for between "
                             f"{SEARCH_LO_FACTOR:g}x and {SEARCH_HI_FACTOR:g}x this value; "
                             "without it the search runs unbounded and can walk into the ambient tail.")
    parser.add_argument("--grid-points", type=int, default=GRID_POINTS, help="Points on the uniform log10(rank) grid the curve is resampled onto.")
    parser.add_argument("--smooth", type=float, default=SMOOTH_FRACTION, help="Moving-average window, as a fraction of the grid.")
    parser.add_argument("--min-umis", type=int, default=100, help="Minimum UMI threshold used to retain entries before analysis.")
    parser.add_argument("--min-points", type=int, default=30, help="Minimum number of retained points required for derivative-based cutoff detection.")
    parser.add_argument("--fallback-cutoff", type=int, default=0, help="Cutoff written when the dataset is too small to analyze.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    try:
        data = np.loadtxt(args.input, dtype=int)
        data = np.atleast_1d(data).flatten()

        umis = np.sort(data)[::-1]
        umis = umis[umis >= args.min_umis]

        n = len(umis)

        if n < 10:
            msg = f"Too few points above {args.min_umis} UMIs for reliable cutoff detection ({n} found)."
            print(f"Warning: {msg}", file=sys.stderr)
            write_fallback(args.sample_id, args.out_json, args.out_cutoff, msg, cutoff=args.fallback_cutoff)
            return

        n_grid = max(50, min(args.grid_points, n))
        smooth_window = max(3, int(n_grid * args.smooth) | 1)

        ranks = np.arange(1, n + 1)
        cum_sum = np.cumsum(umis, dtype=np.float64)
        log_x, log_y = np.log10(ranks), np.log10(cum_sum)

        grid_x, grid_y = resample_log_rank(log_x, log_y, n_grid, smooth_window)
        _, d2 = calculate_derivatives(grid_x, grid_y)

        grid_ranks = 10.0 ** grid_x
        idx, bounded = detect_inflection_point(
            d2, grid_ranks, smooth_window, args.expected_cells,
            SEARCH_LO_FACTOR, SEARCH_HI_FACTOR,
        )

        # Back from the grid to an actual barcode rank, then to the UMI count there.
        cutoff_rank = int(np.clip(round(grid_ranks[idx]), 1, n))
        final_cutoff = int(umis[cutoff_rank - 1])

        status = "ok" if n >= args.min_points else "warning"
        message = (
            "Cutoff calculated successfully."
            if status == "ok"
            else f"Cutoff calculated on small dataset ({n} points above {args.min_umis} UMIs); interpret cautiously."
        )
        if not bounded and args.expected_cells:
            status = "warning"
            message = (
                f"No knee found within {SEARCH_LO_FACTOR:g}-{SEARCH_HI_FACTOR:g}x the "
                f"{args.expected_cells} expected cells; searched the whole curve instead."
            )
        elif not args.expected_cells:
            status = "warning"
            message = (
                "No expected cell count given, so the knee was searched for across the "
                "whole curve; on a sample with a deep ambient tail that can select a "
                "cutoff far below the real cells."
            )

        export_data = {
            args.sample_id: {
                "logX": log_x.tolist(),
                "logY": log_y.tolist(),
                "customdata": [[int(r), int(u)] for r, u in zip(ranks, umis)],
                "derivX": grid_x.tolist(),
                "derivY": (d2 - np.max(d2)).tolist(),
                "threshold_logX": float(log_x[cutoff_rank - 1]),
                "threshold_logY": float(log_y[cutoff_rank - 1]),
                "threshold_umi": final_cutoff,
                "status": status,
                "message": message,
            }
        }

        with open(args.out_json, "w") as f:
            json.dump(export_data, f, separators=(",", ":"))

        with open(args.out_cutoff, "w") as f:
            f.write(str(final_cutoff))

        print(
            f"Processed {args.sample_id} -> UMI Cutoff: {final_cutoff} "
            f"(rank {cutoff_rank} of {n} barcodes above {args.min_umis} UMIs; "
            f"expected_cells={args.expected_cells}, search {'bounded' if bounded else 'unbounded'})",
            file=sys.stderr,
        )

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
