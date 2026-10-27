#!/usr/bin/env python3
"""
Check the added rRNA reference against the base genome (PREPARE_ADDFEATURE).

``ref_fasta_addfeature`` holds extra sequences appended to ``ref_fasta`` for
mapping; ``ref_gtf_addfeature`` annotates rRNA loci, on contigs of either FASTA,
for the read-level rRNA metrics only. Rfam-derived inputs commonly annotate loci
of the genome itself and carry the same loci, cut out of the same assembly, as
the added sequences. Appending such copies turns every read there into a
multimapper, so exact copies of reference sequence are dropped here.

Writes three files:

``--out-fasta``
    The added records that are not copies of reference sequence (may be empty).

``--out-regions``
    The rRNA regions CALC_READ_METRICS counts, as a TSV with a header,
    0-based half-open: ``contig start end name kind``. ``kind`` is
      feature           an added GTF row (any feature type)
      added_contig      the full span of a kept added record
      duplicate_contig  the full span of a dropped copy; still listed so BAMs
                        mapped against a reference that did hold the copy count it

``--out-report``
    One row per added record and per dropped GTF contig: what was done and why.

A record is a copy when its sequence equals a base contig, an earlier added
record, or the locus of the added GTF row whose gene_id or transcript_id is its
name (reverse-complemented on the minus strand). Copies of reference sequence
that no GTF row links by name are not detected.
"""

import argparse
import gzip
import hashlib
import re
import sys
from collections import defaultdict


_ID_RE = re.compile(r'(gene_id|transcript_id) "([^"]*)"')
_COMPLEMENT = str.maketrans("ACGTRYKMBDHVN", "TGCAYRMKVHDBN")

KINDS = ("feature", "added_contig", "duplicate_contig")


def warn(*lines):
    for line in lines:
        print(f"WARNING: {line}", file=sys.stderr)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def digest(seq):
    return hashlib.sha1(seq.encode()).hexdigest()


# ------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------

def read_added_fasta(path):
    """``{name: upper-case sequence}`` in file order; names are the first header word."""
    records, name, chunks = {}, None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                words = line[1:].split()
                name, chunks = (words[0] if words else ""), []
                if not name or name in records:
                    sys.exit(f"ERROR: empty or repeated header '{line}' in {path}")
            elif line:
                if name is None:
                    sys.exit(f"ERROR: sequence before the first header in {path}")
                chunks.append(line.upper())
    if name is not None:
        records[name] = "".join(chunks)
    return records


def read_added_gtf(path):
    """GTF rows as ``(contig, start0, end, strand, name, ids)``.

    *name* is the gene_id, else the transcript_id, else the file position; *ids*
    is the set of both, which links a row to the added record of that name.
    """
    rows = []
    with open_text(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 9:
                sys.exit(f"ERROR: line {n} of {path} has {len(cols)} tab-separated columns, not 9")
            try:
                start, end = int(cols[3]) - 1, int(cols[4])
            except ValueError:
                sys.exit(f"ERROR: line {n} of {path} has non-integer coordinates")
            ids = dict(_ID_RE.findall(cols[8]))
            name = ids.get("gene_id") or ids.get("transcript_id") or f"{path}:{n}"
            rows.append((cols[0], start, end, cols[6], name, set(ids.values())))
    return rows


def scan_base_fasta(path, wanted):
    """One streaming pass over the base genome.

    *wanted* is ``{contig: [(start0, end, key)]}``. Returns ``(lengths, hashes,
    slices)``: every contig's length, the SHA-1 of every contig's upper-case
    sequence, and ``{key: sequence}`` for the wanted intervals. Only those slices
    are held, never a whole chromosome.
    """
    lengths, hashes, slices = {}, {}, defaultdict(list)
    name, pos, sha = None, 0, None
    pending, active, nxt = [], [], 0

    def close():
        if name is not None:
            lengths[name] = pos
            hashes[name] = sha.hexdigest()

    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                close()
                name, pos, sha = line[1:].split()[0], 0, hashlib.sha1()
                pending, active, nxt = sorted(wanted.get(name, ())), [], 0
                continue
            seq = line.strip().upper()
            if not seq:
                continue
            sha.update(seq.encode())
            end = pos + len(seq)
            while nxt < len(pending) and pending[nxt][0] < end:
                active.append(pending[nxt])
                nxt += 1
            if active:
                for istart, iend, key in active:
                    lo, hi = max(istart, pos), min(iend, end)
                    if lo < hi:
                        slices[key].append(seq[lo - pos:hi - pos])
                active = [iv for iv in active if iv[1] > end]
            pos = end
    close()
    return lengths, hashes, {key: "".join(parts) for key, parts in slices.items()}


# ------------------------------------------------------------------
# Classification
# ------------------------------------------------------------------

def linked_loci(records, rows):
    """``{record: [(contig, start0, end, strand)]}``: GTF loci named like an added record."""
    loci = defaultdict(set)
    for contig, start, end, strand, _name, ids in rows:
        for rid in ids & records.keys():
            loci[rid].add((contig, start, end, strand))
    return loci


def classify(records, loci, lengths, hashes, slices):
    """``{record: (status, detail)}``; status is ``kept`` or ``duplicate``.

    Exits on a record named like a base contig but differing from it: the
    merged genome would hold two chromosomes of the same name.
    """
    base_by_hash = {h: c for c, h in hashes.items()}
    seen = {}
    result = {}
    for name, seq in records.items():
        h = digest(seq)
        if name in lengths:
            if hashes[name] != h:
                sys.exit(f"ERROR: added record '{name}' has the name of a ref_fasta contig "
                         "but a different sequence; rename it or remove it")
            result[name] = ("duplicate", f"identical to ref_fasta contig {name}")
        elif h in base_by_hash:
            result[name] = ("duplicate", f"identical to ref_fasta contig {base_by_hash[h]}")
        elif h in seen:
            result[name] = ("duplicate", f"identical to added record {seen[h]}")
        else:
            result[name] = classify_linked(name, seq, loci.get(name, ()), lengths, slices)
        seen.setdefault(h, name)
    return result


def classify_linked(name, seq, loci, lengths, slices):
    """Kept or duplicate, from the base loci the added GTF names like *name*."""
    base_loci = sorted(l for l in loci if l[0] in lengths)
    for contig, start, end, strand in base_loci:
        region = slices.get((contig, start, end), "")
        if strand == "-":
            region = revcomp(region)
        if region == seq:
            return "duplicate", f"identical to {contig}:{start + 1}-{end}({strand}), annotated as {name}"
    if base_loci:
        contig, start, end, strand = base_loci[0]
        warn(f"added record {name} differs from the locus annotated under its name "
             f"({contig}:{start + 1}-{end}); kept, but reads there may still multimap")
        return "kept", f"differs from {contig}:{start + 1}-{end}({strand}), annotated as {name}"
    return "kept", "no copy of it found in ref_fasta"


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base-fasta", required=True, help="ref_fasta")
    p.add_argument("--add-fasta", default="", help="ref_fasta_addfeature")
    p.add_argument("--add-gtf", default="", help="ref_gtf_addfeature")
    p.add_argument("--out-fasta", required=True, help="Output: added records to append to ref_fasta")
    p.add_argument("--out-regions", required=True, help="Output: rRNA regions for CALC_READ_METRICS")
    p.add_argument("--out-report", required=True, help="Output: what was kept and dropped, and why")
    return p.parse_args()


def main():
    args = parse_args()
    records = read_added_fasta(args.add_fasta) if args.add_fasta else {}
    rows = read_added_gtf(args.add_gtf) if args.add_gtf else []

    loci = linked_loci(records, rows)
    wanted = defaultdict(set)
    for locs in loci.values():
        for contig, start, end, _strand in locs:
            wanted[contig].add((start, end, (contig, start, end)))
    lengths, hashes, slices = scan_base_fasta(args.base_fasta, wanted)
    status = classify(records, loci, lengths, hashes, slices)

    with open(args.out_fasta, "w") as fh:
        for name, seq in records.items():
            if status[name][0] == "kept":
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")

    # Feature rows on an added contig are already covered by its full span, and
    # rows on a contig of neither FASTA can never hold a read
    kept_added = {n for n, (s, _) in status.items() if s == "kept"}
    unknown = defaultdict(int)
    with open(args.out_regions, "w") as fh:
        fh.write("contig\tstart\tend\tname\tkind\n")
        for contig, start, end, _strand, name, _ids in rows:
            if contig in lengths or contig in kept_added:
                fh.write(f"{contig}\t{start}\t{end}\t{name}\tfeature\n")
            elif contig not in records:
                unknown[contig] += 1
        for name, seq in records.items():
            # A copy carrying a base contig's own name is that contig, which is no rRNA
            if name in lengths:
                continue
            kind = "added_contig" if status[name][0] == "kept" else "duplicate_contig"
            fh.write(f"{name}\t0\t{len(seq)}\t{name}\t{kind}\n")

    with open(args.out_report, "w") as fh:
        fh.write("item\tname\tstatus\tdetail\n")
        for name, (st, detail) in status.items():
            fh.write(f"fasta_record\t{name}\t{st}\t{detail}\n")
        for contig, n in sorted(unknown.items()):
            fh.write(f"gtf_contig\t{contig}\tdropped\t{n} row(s) on a contig of neither "
                     "ref_fasta nor ref_fasta_addfeature\n")

    duplicates = [n for n, (s, _) in status.items() if s == "duplicate"]
    print(f"Added records: {len(records)}, kept {len(records) - len(duplicates)}, "
          f"dropped as copies of reference sequence {len(duplicates)}", file=sys.stderr)
    if duplicates:
        shown = " ".join(duplicates[:20]) + (" ..." if len(duplicates) > 20 else "")
        warn(f"{len(duplicates)} record(s) of ref_fasta_addfeature copy reference sequence "
             f"and are not appended: {shown}")
    if unknown:
        warn(f"{sum(unknown.values())} row(s) of ref_gtf_addfeature are on contigs of neither "
             "FASTA and are ignored: " + " ".join(sorted(unknown)),
             "check that their contig names match ref_fasta")
    print(f"Added GTF rows: {len(rows)}, of which {sum(unknown.values())} ignored", file=sys.stderr)


if __name__ == "__main__":
    main()
