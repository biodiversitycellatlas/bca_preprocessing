#!/usr/bin/env python3
"""
Split paired FASTQs by group using a whitelist and edit-distance (≤2) barcode matching.

The round-1 barcode is matched like its correction tables do 
(smallest Levenshtein distance, A/C/G/T only), a read is kept
when any of its best-matching wells is in the group, and wells are given in Parse's
notation ('A1:C6' block, 'A1-B6' range in plate order, 'C4' single well).
"""

import argparse, gzip, itertools, os
from collections import Counter
from typing import Dict, Iterable, Tuple, Optional, Set, List
import re

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

PLATE_ROWS = "ABCDEFGH"
PLATE_COLS = 12
WELL_RE = re.compile(r"([A-H])(1[012]|[1-9])$")

def parse_well_range(rng_str):
    """
    Expand a well specification on a 96-well plate:
      'A1:C6' is a block, top-left to bottom-right: A1-A6, B1-B6, C1-C6;
      'A1-B6' is a range in plate order: A1-A12, B1-B6;
      'C4' is a single well; selections are joined by commas, e.g. 'A1-A6,B1:D3,C4'.
    Return the wells in plate order, each once. A block or range that runs backwards
    (e.g. 'B3:A1') selects nothing; a malformed well raises ValueError.
    """
    rng_str = rng_str.replace(" ", "").replace(";", ",").upper()

    def row_col(w):
        m = WELL_RE.match(w)
        if not m:
            raise ValueError(f"Not a well of a 96-well plate (A1-H12): '{w}' in '{rng_str}'")
        return PLATE_ROWS.index(m.group(1)), int(m.group(2)) - 1

    selected = set()
    for part in filter(None, rng_str.split(",")):
        if ":" in part:    # block
            (r0, c0), (r1, c1) = (row_col(w) for w in part.split(":", 1))
            selected |= {r * PLATE_COLS + c for r in range(r0, r1 + 1) for c in range(c0, c1 + 1)}
        elif "-" in part:  # range in plate order
            (r0, c0), (r1, c1) = (row_col(w) for w in part.split("-", 1))
            selected |= set(range(r0 * PLATE_COLS + c0, r1 * PLATE_COLS + c1 + 1))
        else:
            r, c = row_col(part)
            selected.add(r * PLATE_COLS + c)
    return [f"{PLATE_ROWS[i // PLATE_COLS]}{i % PLATE_COLS + 1}" for i in sorted(selected)]

def open_fastq(path: str):
    """Open a FASTQ for reading as text, gzipped or not (detected from the magic bytes)."""
    with open(path, "rb") as probe:
        is_gzipped = probe.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gzipped else open(path, "rt")

def deduce_fq2(fq1: str) -> str:
    """Replace last 'R1' with 'R2' to infer R2 filename."""
    parts = fq1.split("R1")
    if len(parts) < 2: raise ValueError("Cannot deduce --fq2 from --fq1, please provide --fq2.")
    return "R1".join(parts[:-1]) + "R2" + parts[-1]

def open_out(path: str, sample_id: str, sample_name: str, read: str):
    """
    Open gzip text output as <sample_id>_group_<sample_name>_<read>.fastq.gz.

    Level 2, as SCIROCKET_DEMUX writes: gzip's default level 9 compressed real reads
    18x slower for files only a fifth smaller, and took most of the run time.
    """
    os.makedirs(path, exist_ok=True)
    fn = os.path.join(path, f"{sample_id}_group_{sample_name}_{read}.fastq.gz")
    return gzip.open(fn, "wt", compresslevel=2), fn

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

# The distance between the read's barcode slice and a
# whitelist barcode (both 8 nt) is the Levenshtein (edit) distance. Between strings of
# equal length, distance 1 is one substitution, and distance 2 is two substitutions or
# one deletion plus one insertion, e.g. a base lost inside bc1 with the next base shifted
# in. The tables hold A/C/G/T sequences only, so a slice with N (or anything else) never
# matches. Every barcode at the smallest distance is returned, ties included.
BASES = "ACGT"

def edit_distance(a: str, b: str) -> int:
    """Levenshtein distance (substitutions, insertions and deletions each cost 1)."""
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]

def best_matches_all(bc: str, wl_seqs: Iterable[str], max_ed: int = 2) -> Tuple[List[str], Optional[int]]:
    """
    Return (best_seqs, best_distance): all whitelist sequences at the smallest edit distance
    from bc, if that is at most max_ed, else ([], None). A barcode with a character other
    than A/C/G/T has no match. Scans the whole whitelist: the reference for the table.
    """
    bc = bc.upper()
    if not bc or any(c not in BASES for c in bc):
        return [], None
    best_d, best_list = None, []
    for w in wl_seqs:
        d = edit_distance(bc, w)
        if best_d is None or d < best_d:
            best_d, best_list = d, [w]
        elif d == best_d:
            best_list.append(w)
    if best_d is None or best_d > max_ed:
        return [], None
    return best_list, best_d

def _neighbours(w: str, max_ed: int) -> Iterable[Tuple[str, int]]:
    """Yield (sequence, edit distance to w) for the A/C/G/T sequences of w's length within max_ed (<= 2)."""
    for d in range(min(max_ed, 2) + 1):
        for positions in itertools.combinations(range(len(w)), d):
            options = [[c for c in BASES if c != w[p]] for p in positions]
            for subs in itertools.product(*options):
                s = list(w)
                for p, c in zip(positions, subs):
                    s[p] = c
                yield "".join(s), d
    if max_ed >= 2:
        # One deletion plus one insertion. Between equal-length strings the edit distance
        # is 0 if equal, 1 if they differ at one position, and 2 otherwise.
        for i in range(len(w)):
            short = w[:i] + w[i + 1:]
            for j in range(len(short) + 1):
                for c in BASES:
                    s = short[:j] + c + short[j:]
                    mism = sum(x != y for x, y in zip(s, w))
                    yield s, min(mism, 2)

def build_match_table(wl_seqs: Iterable[str], max_ed: int = 2) -> Dict[str, Tuple[int, List[str]]]:
    """
    Return {sequence: (best_distance, best_seqs)} for every A/C/G/T sequence within max_ed
    of a whitelist entry: exactly what best_matches_all returns for it, computed once
    instead of once per read. A sequence not in the table has no match.
    """
    if max_ed > 2:
        raise ValueError("max_ed is at most 2")
    table: Dict[str, Tuple[int, List[str]]] = {}
    for w in dict.fromkeys(s.upper() for s in wl_seqs):
        for key, d in _neighbours(w, max_ed):
            hit = table.get(key)
            if hit is None or d < hit[0]:
                table[key] = (d, [w])
            elif d == hit[0] and hit[1][-1] != w:
                hit[1].append(w)
    return table

def best_matches(bc: str, table: Dict[str, Tuple[int, List[str]]]) -> Tuple[List[str], Optional[int]]:
    """(best_seqs, best_distance) for bc from the match table, or ([], None) without a match."""
    hit = table.get(bc.upper())
    return (hit[1], hit[0]) if hit is not None else ([], None)


def fastq_pairs(fq1: str, fq2: str):
    """Yield paired FASTQ records as (r1_lines, r2_lines) 4-tuples; stops at EOF."""
    with open_fastq(fq1) as f1, open_fastq(fq2) as f2:
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
    table = build_match_table(wl_seqs, max_ed=max_ed)

    for r1, r2 in fastq_pairs(fq1, fq2):
        stats["total"] += 1
        s2 = r2[1].strip().upper()
        if len(s2) < b_end:
            stats["too_short"] += 1
            continue

        bc = s2[b_start:b_end]
        if "N" in bc:
            stats["n_in_bc"] += 1

        best_list, dist = best_matches(bc, table)
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
    ap.add_argument("--fq1", required=True, help="Input R1 FASTQ, plain or gzipped.")
    ap.add_argument("--fq2", help="Input R2 FASTQ, plain or gzipped; inferred from R1 if omitted.")
    ap.add_argument("--whitelist", required=True, help="CSV (bci,sequence,uid,well,stype).")
    ap.add_argument("--group", nargs=2, metavar=("SAMPLE", "WELLS"),
                    help="Sample name and wells, in Parse's notation: A1:C6 block, A1-B6 range in plate order, C4 single; joined by commas.")
    ap.add_argument("--output", required=True, help="Output directory.")
    ap.add_argument("--barcode_start", type=int, default=50, help="0-based start index in R2.")
    ap.add_argument("--barcode_end", type=int, default=58, help="Exclusive end index in R2.")
    ap.add_argument("--max_edit_dist", type=int, default=2, choices=[0, 1, 2], help="Max edit (Levenshtein) distance (≤2).")
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
    print(f"Barcodes with 'N' (never matched): {stats['n_in_bc']}\n")
    print(f"Output R1: {out1}")
    print(f"Output R2: {out2}\n")

if __name__ == "__main__":
    main()
