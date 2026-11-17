#!/usr/bin/env python3
"""
Merge each round-1 well's barcodes into one: rewrite bc1 in the barcode read (R2) to the
well's oligo-dT (T) barcode.

Each Parse round-1 well holds an oligo-dT (T) and a random-hexamer (R) RT primer, and in
some kits targeted (X) primers, each with its own bc1 sequence. split-pipe names cells by
well, so these reads form one cell; a mapper given every bc1 sequence would split each
cell in two or more. bc1 is matched like in parsebio_custom_demux.py and Parse's own
tables (smallest edit distance, at most --max_edit_dist, A/C/G/T only): when all best
matches lie in one well, the 8 bc1 bases become that well's T barcode. Reads tied between
wells, unmatched (including bc1 with N) or too short are written unchanged. Only R2 is
read and written; R1 is not needed.

Also writes parse_bc1_T_whitelist.txt, the T barcode of every well in CSV order, as the
round-1 whitelist for the mappers.
"""

import argparse, csv, gzip, os, sys
from collections import Counter, OrderedDict
from typing import Dict, Tuple

from parsebio_custom_demux import best_matches, build_match_table, open_fastq

WHITELIST_NAME = "parse_bc1_T_whitelist.txt"


def load_wells(path: str) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Read split-pipe's CSV (bci,sequence,uid,well,stype).

    Return ({seq: well}, {seq: stype}, {well: canonical seq}), sequences uppercase. The
    canonical barcode of a well is its T barcode, else its first row in bci order.
    Wells keep the CSV's order of first appearance.
    """
    seq_well: Dict[str, str] = {}
    seq_type: Dict[str, str] = {}
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            seq = (row.get("sequence") or "").strip().upper()
            well = (row.get("well") or "").strip().upper()
            if not seq or not well:
                continue
            bci = row.get("bci", "").strip()
            rows.append((int(bci) if bci.isdigit() else len(rows), seq, well, (row.get("stype") or "").strip()))
    if not rows:
        raise ValueError(f"{path}: no rows with a 'sequence' and a 'well' (expected split-pipe's CSV: bci,sequence,uid,well,stype)")

    canonical: Dict[str, str] = OrderedDict()
    for _, seq, well, stype in rows:
        canonical.setdefault(well, None)
    for _, seq, well, stype in sorted(rows):
        seq_well[seq] = well
        seq_type[seq] = stype
        if stype == "T" and canonical[well] is None:
            canonical[well] = seq
    for _, seq, well, _ in sorted(rows):
        if canonical[well] is None:
            canonical[well] = seq
    return seq_well, seq_type, canonical


def collapse(fq2: str, out, seq_well: Dict[str, str], seq_type: Dict[str, str],
             canonical: Dict[str, str], b_start: int, b_end: int, max_ed: int = 2) -> Dict[str, int]:
    """Stream R2 records to out, with bc1 rewritten to its well's canonical barcode."""
    stats = Counter(total=0, rewritten=0, ambiguous=0, unmatched=0, too_short=0)
    wl_seqs = list(seq_well)
    table = build_match_table(wl_seqs, max_ed=max_ed)

    with open_fastq(fq2) as fh:
        while True:
            rec = [fh.readline() for _ in range(4)]
            if not rec[1]:
                break
            stats["total"] += 1
            seq = rec[1].rstrip("\r\n")
            if len(seq) < b_end:
                stats["too_short"] += 1
            else:
                best, dist = best_matches(seq[b_start:b_end], table)
                wells = {seq_well[b] for b in best}
                if not best:
                    stats["unmatched"] += 1
                elif len(wells) > 1:
                    stats["ambiguous"] += 1
                else:
                    stats["rewritten"] += 1
                    stats[f"ed{dist}"] += 1
                    # Primer type(s) of the best matches: from_R, from_T, from_X..., or
                    # e.g. from_R+T for a read tied between two barcodes of one well.
                    types = sorted({seq_type[b] or "?" for b in best})
                    stats[f"from_{'+'.join(types)}"] += 1
                    rec[1] = seq[:b_start] + canonical[wells.pop()] + seq[b_end:] + "\n"
            out.write("".join(rec))
    return dict(stats)


def main():
    ap = argparse.ArgumentParser(description="Rewrite Parse bc1 to each well's T barcode.")
    ap.add_argument("--csv", required=True, help="split-pipe's round-1 CSV (bci,sequence,uid,well,stype).")
    ap.add_argument("--fq2", required=True, help="Barcode read (R2) FASTQ, plain or gzipped.")
    ap.add_argument("--sample_id", required=True, help="Sample ID used in output filenames.")
    ap.add_argument("--output", default=".", help="Output directory.")
    ap.add_argument("--barcode_start", type=int, default=50, help="0-based start of bc1 in R2.")
    ap.add_argument("--barcode_end", type=int, default=58, help="Exclusive end of bc1 in R2.")
    ap.add_argument("--max_edit_dist", type=int, default=2, choices=[0, 1, 2], help="Max edit (Levenshtein) distance (<=2).")
    args = ap.parse_args()

    seq_well, seq_type, canonical = load_wells(args.csv)
    lens = {len(s) for s in seq_well}
    if lens != {args.barcode_end - args.barcode_start}:
        sys.exit(f"bc1 slice length {args.barcode_end - args.barcode_start} does not match the CSV barcodes ({sorted(lens)}).")

    os.makedirs(args.output, exist_ok=True)
    with open(os.path.join(args.output, WHITELIST_NAME), "w") as fh:
        fh.write("\n".join(canonical.values()) + "\n")

    out_name = os.path.join(args.output, f"{args.sample_id}_R2.bc1collapsed.fastq.gz")
    with gzip.open(out_name, "wt", compresslevel=2) as out:
        stats = collapse(args.fq2, out, seq_well, seq_type, canonical,
                         args.barcode_start, args.barcode_end, max_ed=args.max_edit_dist)

    lines = [f"Sample:\t{args.sample_id}",
             f"Wells:\t{len(canonical)} ({len(seq_well)} round-1 barcodes)"]
    lines += [f"{k}:\t{stats[k]}" for k in ("total", "rewritten", "ambiguous", "unmatched", "too_short")]
    lines += [f"{k}:\t{v}" for k, v in sorted(stats.items()) if k.startswith(("ed", "from_"))]
    report = "\n".join(lines) + "\n"
    sys.stdout.write(report)
    with open(os.path.join(args.output, f"{args.sample_id}_bc1_collapse.log"), "w") as fh:
        fh.write(report)


if __name__ == "__main__":
    main()
