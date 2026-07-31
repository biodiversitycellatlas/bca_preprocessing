#!/usr/bin/env python3
"""
Pre-calculate a second-derivative-based cell-calling cutoff.
"""

import argparse
import json
import sys
from typing import Tuple

import numpy as np


def moving_average(data: np.ndarray, window_size: int) -> np.ndarray:
    """
    Smooth the array using a centered moving average.
    """
    if window_size < 1:
        return data
    return np.convolve(data, np.ones(window_size) / window_size, mode="same")


def calculate_derivatives(log_x: np.ndarray, log_y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute first and second discrete derivatives.
    """
    d1 = np.diff(log_y) / np.diff(log_x)
    d2 = np.diff(d1) / np.diff(log_x[:-1])
    return d1, d2


def detect_inflection_point(smoothed_2nd: np.ndarray, trim_edge: int = 15) -> int:
    """
    Detect the index of the most negative second-derivative point.
    Edge trimming avoids unstable extrema near the boundaries.
    """
    if len(smoothed_2nd) <= (trim_edge * 2):
        return 0
    return int(np.argmin(smoothed_2nd[trim_edge:-trim_edge]) + trim_edge)


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
    parser.add_argument("--smooth", type=int, default=200, help="Moving-average window size applied to the second derivative.")
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

        effective_smooth = min(args.smooth, max(3, n // 4))
        if effective_smooth % 2 == 0:
            effective_smooth -= 1
        effective_smooth = max(effective_smooth, 3)

        status = "ok" if n >= args.min_points else "warning"
        message = (
            "Cutoff calculated successfully."
            if status == "ok"
            else f"Cutoff calculated on small dataset ({n} points above {args.min_umis} UMIs); interpret cautiously."
        )

        ranks = np.arange(1, n + 1)
        cum_sum = np.cumsum(umis)
        log_x, log_y = np.log10(ranks), np.log10(cum_sum)

        _, d2 = calculate_derivatives(log_x, log_y)
        smoothed_d2 = moving_average(d2, effective_smooth)

        idx = detect_inflection_point(smoothed_d2)
        final_cutoff = int(umis[idx])

        export_data = {
            args.sample_id: {
                "logX": log_x.tolist(),
                "logY": log_y.tolist(),
                "customdata": [[int(r), int(u)] for r, u in zip(ranks, umis)],
                "derivX": log_x[1:-1].tolist(),
                "derivY": (smoothed_d2 - np.max(smoothed_d2)).tolist(),
                "threshold_logX": float(log_x[idx]),
                "threshold_logY": float(log_y[idx]),
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
            f"Processed {args.sample_id} -> UMI Cutoff: {final_cutoff}",
            file=sys.stderr,
        )

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
