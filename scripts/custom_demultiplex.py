#!/usr/bin/env python3

"""
demux_assign_groups.py

Demultiplex paired-end FASTQs by barcodes (CB+UMI on R2).
We:
  - Parse a 5-column barcode CSV: (bci,sequence,uid,well,stype)
  - Map each (barcode_seq) -> (well)
  - Define user "groups", e.g.: [("group_name", "A1-A2"), ...]
  - Read R1 (cDNA) / R2 (CB+UMI) using umi_tools (paired-end)
  - Extract barcodes from R2 with optional ± frameshift
  - Find matching well from the barcode
  - Find group for that well
  - Write read pair to group-based FASTQs
  - At the end, report how many reads per group, plus how many unassigned.

Usage:
    python demux_assign_groups.py \
        --r1 R1.fastq.gz \
        --r2 R2.fastq.gz \
        --bc-csv barcodes.csv \
        --bc-start 0 \
        --bc-len 8 \
        --frameshift 1 \
        --outdir demux_out

Example groups are hardcoded below (in `GROUP_DEFS`) as you mentioned:
    GROUP_DEFS = [
        ["1010_ACMEsorb_GM",   "A1-A2"],
        ["1010_ACMEsorb_Lib",  "A3-A4"],
        ...
    ]
Adapt as needed.
"""

import sys
import os
import gzip
import argparse
import re
from collections import defaultdict

# umi_tools for FASTQ parsing
from umi_tools import fastqtransform


############################################################################
# 1) Hardcode your group definitions or read from a file
############################################################################
GROUP_DEFS = [
    ["1010_ACMEsorb_GM",   "A1-A2"],
    ["1010_ACMEsorb_Lib",  "A3-A4"],
    ["1010_DSP_CMFSV_Lib", "A5-A6"],
    ["1710_ACMEsorb_GM",   "A7-A8"],
    ["1710_ACMEsorb_Lib",  "A9-A10"],
    ["1710_DSP_CMFSW_Lib", "A11-A12"]
]

def parse_well_range(well_str):
    """
    Minimal parser for well spec like "A1-A2" or single "A1".
    If needed for an 8x12 plate, adapt further.
    E.g., "A1-A6" => ['A1','A2','A3','A4','A5','A6'].
    """
    well_str = well_str.strip().upper()

    # If it's a range "A1-A2"
    if "-" in well_str:
        start, end = well_str.split("-")
        row_s, col_s = re.match(r"([A-H])(\d+)", start).groups()
        row_e, col_e = re.match(r"([A-H])(\d+)", end).groups()
        col_s, col_e = int(col_s), int(col_e)
        # For simplicity, assume row_s == row_e
        out = []
        step = 1 if col_s <= col_e else -1
        for c in range(col_s, col_e + step, step):
            out.append(f"{row_s}{c}")
        return out
    else:
        # single well
        return [well_str]

# Build a global dictionary well->group
WELL2GROUP = {}
for gname, wspec in GROUP_DEFS:
    wells = parse_well_range(wspec)
    for w in wells:
        WELL2GROUP[w] = gname


############################################################################
# 2) Read the 5-column barcode CSV, store as { barcode_seq -> well }
############################################################################
def load_barcode_map(bc_csv):
    """
    Expects columns: bci,sequence,uid,well,stype
    We'll store bc_map[sequence] = well.
    If you need to filter by stype, do so here.
    """
    bc_map = {}
    with open(bc_csv, "r") as f:
        # Skip a header if it exists. Check line 1:
        first_line = f.readline().strip()
        # If it doesn't match expected columns, treat as data
        if not first_line.lower().startswith("bci,sequence,uid,well,stype"):
            # It's data, so parse it
            row = first_line.split(",")
            if len(row) == 5:
                # parse
                bci, seq, uid, well, stype = row
                bc_map[seq.strip()] = well.strip()
            # else might be malformed
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) < 5:
                continue
            bci, seq, uid, well, stype = parts[:5]
            seq  = seq.strip()
            well = well.strip()
            bc_map[seq] = well
    return bc_map


############################################################################
# 3) Barcode extraction with frameshift
############################################################################
def get_bc_candidates(seq, start, length, frameshift=0):
    """
    Return all possible substrings of length 'length'
    from 'seq' at positions [start - frameshift .. start + frameshift].
    """
    cands = []
    for shift in range(-frameshift, frameshift+1):
        s0 = start + shift
        e0 = s0 + length
        if s0 < 0:
            continue
        if e0 > len(seq):
            continue
        cands.append(seq[s0:e0])
    return cands


############################################################################
# Main CLI
############################################################################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True, help="R1 FASTQ.gz (cDNA)")
    parser.add_argument("--r2", required=True, help="R2 FASTQ.gz (CB+UMI)")
    parser.add_argument("--bc-csv", required=True, help="5-column barcode CSV (bci,sequence,uid,well,stype)")
    parser.add_argument("--bc-start", type=int, default=0, help="Barcode start on R2")
    parser.add_argument("--bc-len",   type=int, default=8, help="Barcode length on R2")
    parser.add_argument("--frameshift", type=int, default=0, help="Allow ± this many bases shift around bc-start")
    parser.add_argument("--outdir", default="demux_out", help="Output directory")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1) Load barcodes => { sequence -> well }
    bc_map = load_barcode_map(args.bc_csv)
    print(f"Loaded {len(bc_map)} barcodes from {args.bc_csv}")

    # 2) We'll open one pair of output files per group name
    #    plus track counts in assigned_counts
    group_handles = {}
    assigned_counts = defaultdict(int)
    unassigned_count = 0
    total = 0

    def get_handles_for_group(gname):
        if gname not in group_handles:
            r1_out = os.path.join(args.outdir, f"{gname}_R1.fastq.gz")
            r2_out = os.path.join(args.outdir, f"{gname}_R2.fastq.gz")
            fh1 = gzip.open(r1_out, "wb")
            fh2 = gzip.open(r2_out, "wb")
            group_handles[gname] = (fh1, fh2)
        return group_handles[gname]

    # 3) Iterate over paired-end reads using umi_tools function
    pair_iter = fastqtransform._fastq_iter2(args.r1, args.r2)

    for read1, read2 in pair_iter:
        total += 1
        r2seq = read2.seq
        # gather candidate barcodes for frameshift
        bc_cands = get_bc_candidates(r2seq, args.bc_start, args.bc_len, frameshift=args.frameshift)

        matched_well = None
        for c in bc_cands:
            if c in bc_map:
                matched_well = bc_map[c]
                break

        if not matched_well:
            unassigned_count += 1
            continue

        # Which group is that well in?
        gname = WELL2GROUP.get(matched_well, None)
        if not gname:
            unassigned_count += 1
            continue

        # Write the read pair to that group's FASTQs
        fh1, fh2 = get_handles_for_group(gname)

        # R1
        fh1.write(f"@{read1.name}\n".encode())
        fh1.write(f"{read1.seq}\n".encode())
        fh1.write("+\n".encode())
        fh1.write(f"{read1.qual}\n".encode())

        # R2
        fh2.write(f"@{read2.name}\n".encode())
        fh2.write(f"{read2.seq}\n".encode())
        fh2.write("+\n".encode())
        fh2.write(f"{read2.qual}\n".encode())

        assigned_counts[gname] += 1

    # 4) Close files
    for (fh1, fh2) in group_handles.values():
        fh1.close()
        fh2.close()

    # 5) Report stats
    print(f"\nTotal read pairs processed: {total}")
    for g in sorted(assigned_counts.keys()):
        print(f"Group {g} => {assigned_counts[g]} read pairs")
    print(f"Unassigned => {unassigned_count} read pairs\n")


if __name__ == "__main__":
    main()
