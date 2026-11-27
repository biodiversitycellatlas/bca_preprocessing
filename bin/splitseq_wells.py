#!/usr/bin/env python3
"""
SPLiT-seq barcode lists for the mappers, and the round-1 well CSV.

'bc_whitelist' gives one list for all three barcode rounds (the 96 SPLiT-seq barcodes of
Rosenberg et al. 2018), or three lists in the order round 1 (RT), round 2, round 3. A list
may hold bare barcodes or delimited rows, e.g. split-seq-pipeline's 'index,sequence': the
one field made only of A/C/G/T is the barcode. A first row without such a field is a
header, lines starting with '#' are skipped, and gzipped lists are read as well.

Writes:
  splitseq_bc1.txt, splitseq_bc2.txt, splitseq_bc3.txt   one barcode per line
  splitseq_wells.csv                                     bci,sequence,uid,well,stype

The well CSV is read by parsebio_custom_demux.py and parsebio_collapse_bc1.py. Unless
--wells gives one, it is derived from the round-1 list: each RT well holds oligo-dT primer
i (1-48) and random-hexamer primer i+48, so barcodes i and i+48 share well i, numbered
A1-D12 in plate order, with stype T and R.
"""

import argparse, csv, gzip, os, re, shutil, sys
from typing import List

SPLIT = re.compile(r"[,;\t ]+")
DNA = re.compile(r"[ACGT]+")


def read_list(path: str) -> List[str]:
    """Return the barcodes of a list, uppercase, in file order; exit on a malformed list."""
    with open(path, "rb") as probe:
        gzipped = probe.read(2) == b"\x1f\x8b"
    seqs: List[str] = []
    rows = 0
    with (gzip.open(path, "rt") if gzipped else open(path)) as fh:
        for n, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = [f for f in SPLIT.split(line.upper()) if DNA.fullmatch(f)]
            rows += 1
            if not fields and rows == 1:
                continue                    # header
            if len(fields) != 1:
                sys.exit(f"{path}, line {n}: expected one barcode of A/C/G/T, found {len(fields)}: '{line}'")
            seqs.append(fields[0])
    if not seqs:
        sys.exit(f"{path}: no barcodes")
    lens = {len(s) for s in seqs}
    if len(lens) != 1:
        sys.exit(f"{path}: barcodes of different lengths {sorted(lens)}")
    seen, dups = set(), []
    for s in seqs:
        if s in seen:
            dups.append(s)
        seen.add(s)
    if dups:
        sys.exit(f"{path}: barcode(s) listed twice: {', '.join(dups[:5])}")
    return seqs


def derive_wells(bc1: List[str], source: str, out_path: str) -> None:
    """Write the well CSV of the 96 SPLiT-seq RT barcodes."""
    if len(bc1) != 96:
        sys.exit(f"{source}: {len(bc1)} round-1 barcodes. The wells are derived from the 96 SPLiT-seq RT barcodes "
                 "(1-48 oligo-dT, 49-96 random hexamer); for another plate layout set 'splitseq_wells'.")
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["bci", "sequence", "uid", "well", "stype"])
        for i, seq in enumerate(bc1):
            k = i % 48
            well = f"{'ABCD'[k // 12]}{k % 12 + 1}"
            stype = "T" if i < 48 else "R"
            w.writerow([i + 1, seq, f"{well}_{stype}", well, stype])


def check_wells(path: str) -> int:
    """Check a given well CSV the way the demultiplexing reads it; return its row count."""
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh)
                if (r.get("sequence") or "").strip() and (r.get("well") or "").strip()]
    if not rows:
        sys.exit(f"'splitseq_wells' ({path}): no rows with a 'sequence' and a 'well' (expected bci,sequence,uid,well,stype)")
    return len(rows)


def main():
    ap = argparse.ArgumentParser(description="SPLiT-seq barcode lists and round-1 well CSV.")
    ap.add_argument("lists", nargs="+", help="One list for all rounds, or the round 1, 2 and 3 lists.")
    ap.add_argument("--wells", help="Round-1 well CSV to use instead of deriving it.")
    ap.add_argument("--output", default=".", help="Output directory.")
    args = ap.parse_args()

    if len(args.lists) not in (1, 3):
        sys.exit(f"'bc_whitelist' gives {len(args.lists)} files: give one list for all three barcode rounds, "
                 "or three lists in the order round 1 (RT), round 2, round 3.")
    segments = [read_list(p) for p in args.lists]
    if len(segments) == 1:
        segments *= 3

    os.makedirs(args.output, exist_ok=True)
    for k, seqs in enumerate(segments, 1):
        with open(os.path.join(args.output, f"splitseq_bc{k}.txt"), "w") as out:
            out.write("\n".join(seqs) + "\n")

    wells_out = os.path.join(args.output, "splitseq_wells.csv")
    if args.wells:
        n = check_wells(args.wells)
        if os.path.abspath(args.wells) != os.path.abspath(wells_out):
            shutil.copyfile(args.wells, wells_out)
        print(f"Round-1 wells from 'splitseq_wells': {n} barcodes")
    else:
        derive_wells(segments[0], args.lists[0], wells_out)
        print("Round-1 wells derived from the 96 RT barcodes: A1-D12, oligo-dT 1-48, random hexamer 49-96")
    for k, seqs in enumerate(segments, 1):
        print(f"Round {k}: {len(seqs)} barcodes of {len(seqs[0])} nt")


if __name__ == "__main__":
    main()
