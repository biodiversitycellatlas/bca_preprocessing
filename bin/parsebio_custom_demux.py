#!/usr/bin/env python3
"""
Split paired FASTQs by group using a whitelist and edit-distance (≤2) barcode matching.
"""

import argparse, gzip, os
from collections import Counter
from typing import Dict, Iterable, Tuple, Optional, Set, List

# --------------------------- Utilities ---------------------------

def load_whitelist_csv(path: str) -> Dict[str, str]:
    """Load CSV (bci,sequence,uid,well,stype) → dict {SEQ: WELL} in uppercase."""
    wl: Dict[str, str] = {}
    with open(path, "r") as f:
        _ = f.readline()
        for line in f:
            p = line.strip().split(",")
            if len(p) >= 5 and p[1] and p[3]:
                wl[p[1].strip().upper()] = p[3].strip().upper()
    return wl

def parse_well_range(spec: str) -> Set[str]:
    """Parse same-row range 'A1-A3' or single well 'B6' → set of wells."""
    spec = spec.strip().upper()
    if "-" in spec:
        s, e = spec.split("-")
        if s[0] != e[0]: raise ValueError("Only same-row ranges supported (e.g., A1-A6).")
        a, b = int(s[1:]), int(e[1:])
        return {f"{s[0]}{i}" for i in range(a, b + 1)}
    return {spec}

def deduce_fq2(fq1: str) -> str:
    """Replace last 'R1' with 'R2' to infer R2 filename."""
    parts = fq1.split("R1")
    if len(parts) < 2: raise ValueError("Cannot deduce --fq2 from --fq1, please provide --fq2.")
    return "R1".join(parts[:-1]) + "R2" + parts[-1]

def open_out(path: str, sample_id: str, sample_name: str, read: str):
    """Open gzip text output as <sample_id>_group_<sample_name>_<read>.fastq.gz."""
    os.makedirs(path, exist_ok=True)
    fn = os.path.join(path, f"{sample_id}_group_{sample_name}_{read}.fastq.gz")
    return gzip.open(fn, "wt"), fn

# --------------------- Chemistry → bc1 bounds ---------------------

def amp_seq_for_chem(chem: str) -> str:
    """Return amplicon template string for v1/v2 vs v3 (digits mark bc rounds)."""
    c = chem.lower().strip()
    if c in ("v3", "3"):
        return "NNNNNNNNNN33333333ATGAGGGGTCAG22222222TCCAACCACCTC11111111"
    # v1/v2 default
    return "NNNNNNNNNN33333333GTGGCCGATGTTTCGCATCGGCGTACGACT22222222ATCCACGTGCTTGAGACTGTGG11111111"

def part_bounds(template: str) -> Tuple[List[int], List[int]]:
    """Return (starts, ends) for parts index: [0]=polyN, [1]=bc1, [2]=bc2, [3]=bc3]."""
    DNA = set("ACGT")
    non_dna = set(template) - DNA
    parts = len(non_dna) + (0 if "N" in non_dna else 1)
    starts, ends, prev = [-1]*parts, [-1]*parts, -1
    for i, ch in enumerate(template):
        idx = (-1 if ch in DNA else (0 if ch == "N" else int(ch)))
        if idx != prev:
            if idx >= 0: starts[idx] = i
            if prev >= 0: ends[prev] = i
        prev = idx
    if prev >= 0 and ends[prev] < 0: ends[prev] = len(template)
    return starts, ends

def bc1_bounds_from_chem(chem: Optional[str]) -> Optional[Tuple[int, int]]:
    """Return (start,end) for bc1 if chemistry provided, else None."""
    if not chem: return None
    t = amp_seq_for_chem(chem)
    s, e = part_bounds(t)
    return s[1], e[1]

# ------------------------ Edit-distance matching ----------------------------

def hamming(a: str, b: str) -> int:
    """Return Hamming distance for equal-length strings, large int if lengths differ."""
    if len(a) != len(b): return 10**9
    return sum(x != y for x, y in zip(a, b))

def best_matches_all(bc: str, wl_seqs: Iterable[str], max_ed: int = 2) -> Tuple[List[str], Optional[int]]:
    """
    Return (list_of_best_seqs, best_distance) where list_of_best_seqs are
    all whitelist sequences at the minimal Hamming distance ≤ max_ed.
    If no match within max_ed, returns ([], None).
    """
    bc = bc.upper()
    best_d = None
    best_list: List[str] = []
    for w in wl_seqs:
        d = hamming(bc, w)
        if best_d is None or d < best_d:
            best_d = d
            best_list = [w]
        elif d == best_d:
            best_list.append(w)
    if best_d is None or best_d > max_ed:
        return [], None
    return best_list, best_d


def fastq_pairs(fq1: str, fq2: str):
    """Yield paired FASTQ records as (r1_lines, r2_lines) 4-tuples; stops at EOF."""
    with gzip.open(fq1, "rt") as f1, gzip.open(fq2, "rt") as f2:
        while True:
            r1 = [f1.readline() for _ in range(4)]
            r2 = [f2.readline() for _ in range(4)]
            if not r1[1] or not r2[1]: return
            yield r1, r2

# --------------------------- Core ---------------------------

def process_fastq(
    fq1: str,
    fq2: str,
    wl_full: Dict[str, str],        # {barcode_seq: well_str}
    group_wells: Optional[Set[str]], # wells to keep, e.g. {"A1","A2"}; None = all wells
    b_start: int,
    b_end: int,
    out_r1,
    out_r2,
    max_ed: int = 2,
) -> Dict[str, int]:
    """
    Stream pairs, slice barcode from R2, match against *full* whitelist (≤max_ed),
    and write reads if their best barcodes intersect group_wells.
    """
    stats = Counter(
        total=0,
        too_short=0,
        matched=0,
        unmatched=0,
        ed0=0,
        ed1=0,
        ed2=0,
        n_in_bc=0,
        no_group_match=0,   # valid barcode, but no best match in this group
    )
    wl_seqs = list(wl_full.keys())

    for r1, r2 in fastq_pairs(fq1, fq2):
        stats["total"] += 1
        s2 = r2[1].strip().upper()
        if len(s2) < b_end:
            stats["too_short"] += 1
            continue

        bc = s2[b_start:b_end]
        if "N" in bc:
            stats["n_in_bc"] += 1

        best_list, dist = best_matches_all(bc, wl_seqs, max_ed=max_ed)
        if not best_list:
            stats["unmatched"] += 1
            continue

        # Map best barcodes to wells
        wells = {wl_full[b] for b in best_list}
        # Does this read belong to this group
        if group_wells is not None and wells.isdisjoint(group_wells):
            stats["no_group_match"] += 1
            continue

        stats["matched"] += 1
        if dist in (0, 1, 2):
            stats[f"ed{dist}"] += 1

        out_r1.write("".join(r1))
        out_r2.write("".join(r2))

    return dict(stats)


# --------------------------- CLI ---------------------------

def main():
    """Parse args, run splitter, print stats and output paths."""
    ap = argparse.ArgumentParser(description="Split FASTQs by group with edit-distance barcode matching (≤2).")
    ap.add_argument("--sample_id", required=True, help="Sample ID used in output filenames.")
    ap.add_argument("--fq1", required=True, help="Input R1 FASTQ.GZ.")
    ap.add_argument("--fq2", help="Input R2 FASTQ.GZ; inferred from R1 if omitted.")
    ap.add_argument("--whitelist", required=True, help="CSV (bci,sequence,uid,well,stype).")
    ap.add_argument("--group", nargs=2, metavar=("SAMPLE", "WELLS"),
                    help="Sample name and well range (same row), e.g., MYRUN A1-A6.")
    ap.add_argument("--output", required=True, help="Output directory.")
    ap.add_argument("--barcode_start", type=int, default=49, help="0-based start index in R2.")
    ap.add_argument("--barcode_end", type=int, default=57, help="Exclusive end index in R2.")
    ap.add_argument("--max_edit_dist", type=int, default=2, choices=[0, 1, 2], help="Max Hamming distance (≤2).")
    ap.add_argument("--chemistry", choices=["v1", "v2", "v3"], help="Optional: auto-derive bc1 bounds for chemistry.")
    args = ap.parse_args()

    wl_full = load_whitelist_csv(args.whitelist)

    if args.group:
        sample_name, wells = args.group
        keep = parse_well_range(wells)   # e.g. {"A1", "A2", ...}
        group_wells = keep
        if not group_wells:
            raise ValueError(f"No wells parsed from group spec: {wells}")
    else:
        sample_name = "all_samples"
        group_wells = None  # meaning "all wells"

    # Compute bc1 bounds from chemistry if provided; else use explicit indices.
    if args.chemistry:
        b_s, b_e = bc1_bounds_from_chem(args.chemistry)
    else:
        b_s, b_e = args.barcode_start, args.barcode_end

    # Sanity: whitelist barcode length vs chosen slice length.
    bc_lens = {len(s) for s in wl_full.keys()}
    if len(bc_lens) == 1 and (b_e - b_s) != list(bc_lens)[0]:
        print(f"# WARNING: slice length {(b_e - b_s)} ≠ whitelist length {list(bc_lens)[0]} (check chemistry/indices).")

    fq2 = args.fq2 or deduce_fq2(args.fq1)
    out_r1, out1 = open_out(args.output, args.sample_id, sample_name, "R1")
    out_r2, out2 = open_out(args.output, args.sample_id, sample_name, "R2")
    try:
        stats = process_fastq(
            args.fq1,
            fq2,
            wl_full=wl_full,
            group_wells=group_wells,
            b_start=b_s,
            b_end=b_e,
            out_r1=out_r1,
            out_r2=out_r2,
            max_ed=args.max_edit_dist,
        )
    finally:
        out_r1.close()
        out_r2.close()


    print("\nProcessing complete.\n")
    print(f"Total reads:                 {stats['total']}")
    print(f"Too short (R2 < end):        {stats['too_short']}")
    print(f"Matched (≤{args.max_edit_dist}):          {stats['matched']}")
    print(f"  - edit distance 0:         {stats.get('ed0', 0)}")
    print(f"  - edit distance 1:         {stats.get('ed1', 0)}")
    print(f"  - edit distance 2:         {stats.get('ed2', 0)}")
    print(f"Unmatched:                   {stats['unmatched']}")
    print(f"Barcodes with 'N':           {stats['n_in_bc']}\n")
    print(f"Output R1: {out1}")
    print(f"Output R2: {out2}\n")

if __name__ == "__main__":
    main()
