#!/usr/bin/env python3
"""
Barcode demultiplexing for ScaleBio single-cell RNA kits (QuantumScale RNA and
ScaleRNA v1.1 / Extended Throughput).

Three subcommands:

  prepare  Read the user's barcode lists (one per barcode segment, in the layout's
           segment order), check them, and write
             - barcodes.json: every segment's entries (id, sequences, canonical sequence);
             - whitelist_<n>_<segment>.txt: the canonical sequences, one per line, which
               the mappers use as the per-segment whitelists of the synthetic barcode read.
  demux    Stream one chunk of a library (one lane: cDNA + barcode read + index reads; or,
           with --cram, an Ultima CRAM whose barcodes the Ultima trimmer moved into SAM
           tags), correct every barcode against its list (Hamming distance), trim the cDNA,
           assign the read to a sample by its RT well (and PCR well) and write per sample
             <sample>.<chunk>_cDNA.fastq.gz    the trimmed cDNA read
             <sample>.<chunk>_BC_UMI.fastq.gz  a synthetic read: the cell barcode, then the
                                               UMI. The cell barcode is the corrected
                                               canonical barcodes in segment order (v1.1),
                                               or, for QuantumScale, the cell's combination
                                               of list entries written in 16 nt
           plus <chunk>.metrics.json.
  gather   Sum the chunk metrics of the libraries a sample was read from into one
           JSON for that sample.
  decode   Cell barcodes back to their barcode ids (wells, beads) per segment.

Barcode list format: headerless, tab- or whitespace-separated 'id seq [seq ...]'; an
entry with more than one sequence (e.g. a PCR well holding a pool of i5 indices) is one
barcode whose first sequence is its canonical one. A one-column file of sequences is
also read, with the 1-based line number as id.
"""

import argparse
import gzip
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter

# ---------------------------------------------------------------------------
# Read layouts
# ---------------------------------------------------------------------------
#
# Each barcode segment is read from one of the reads 'bc' (the barcode/UMI read: R2 for
# QuantumScale, R1 for v1.1), 'i1' or 'i2', at a fixed position or, with 'anchor', at an
# offset from the splint linker found in the barcode read. 'width' is the segment's width
# in the synthetic read (shorter canonical sequences are padded to it). The synthetic read
# is the segments in 'segments' order, then the UMI; conf/seqtech_parameters.config gives
# STARsolo and alevin-fry the same positions.
LAYOUTS = {
    "quantum": {
        "segments": ["bead1", "bead2", "bead3", "rt", "pcr"],
        "extract": {
            "bead1": ("i1", 0, 8),
            "bead2": ("i1", 12, 8),
            "bead3": ("i1", 24, 8),
            "rt":    ("bc", 8, 8),
            "pcr":   ("i2", 0, 8),
        },
        "width": {"bead1": 8, "bead2": 8, "bead3": 8, "rt": 8, "pcr": 8},
        # The 40 nt of barcodes exceed STARsolo's 31-nt cell barcode, and the ~10^9
        # combinations its whitelists, so the cell is written as its combination of list
        # entries in 16 base-4 digits instead (see encode_cell; 'decode' reverses it)
        "cb_code": 16,
        "umi": ("bc", 0, 8),
        "linker": None,
        "index_reads": ["i1", "i2"],
        "trim": {"A": 8, "G": 12},
    },
    "v1.1": {
        "segments": ["lig", "rt", "pcr"],
        # lig is the barcode read up to the splint; umi and rt follow the splint
        "extract": {
            "lig": ("bc", 0, None),
            "rt":  ("anchor", 15, 10),
            "pcr": ("i2", 0, 10),
        },
        "width": {"lig": 9, "rt": 10, "pcr": 10},
        "umi": ("anchor", 7, 8),
        "linker": {"sequence": "TCAGAGC", "starts": [8, 9], "max_mismatches": 1},
        "index_reads": ["i2"],
        "trim": {"A": 8, "T": 8},
    },
}

PAD_BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def cb_length(layout):
    lay = LAYOUTS[layout]
    return lay.get("cb_code") or sum(lay["width"].values())


def synthetic_length(layout):
    return cb_length(layout) + LAYOUTS[layout]["umi"][2]


def encode_cell(indices, sizes, digits):
    """
    A cell's barcode entries (one index per segment, in segment order) as one number in
    mixed radix (the segments' list sizes), written in 'digits' base-4 digits A/C/G/T,
    most significant first.
    """
    n = 0
    for i, size in zip(indices, sizes):
        n = n * size + i
    out = []
    for _ in range(digits):
        n, d = divmod(n, 4)
        out.append("ACGT"[d])
    return "".join(reversed(out))


def decode_cell(code, sizes):
    n = 0
    for c in code:
        n = n * 4 + "ACGT".index(c)
    indices = []
    for size in reversed(sizes):
        n, i = divmod(n, size)
        indices.append(i)
    if n:
        raise ValueError(f"{code} is beyond the barcode lists")
    return list(reversed(indices))


# ---------------------------------------------------------------------------
# Barcode lists
# ---------------------------------------------------------------------------

def read_barcode_list(path):
    """[(id, [seq, ...]), ...] in file order."""
    entries = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for lineno, line in enumerate(fh, 1):
            fields = line.strip().split()
            if not fields or fields[0].startswith("#"):
                continue
            if len(fields) == 1:
                bc_id, seqs = str(lineno), fields
            else:
                bc_id, seqs = fields[0], fields[1:]
            seqs = [s.upper() for s in seqs]
            bad = [s for s in seqs if not re.fullmatch(r"[ACGT]+", s)]
            if bad:
                raise ValueError(f"{path}:{lineno}: not a barcode sequence: {bad[0]!r}")
            entries.append((bc_id, seqs))
    if not entries:
        raise ValueError(f"{path}: no barcodes")
    return entries


def prepare_segment(name, entries, width, path):
    """Check one segment's list and choose its canonical sequences, padded to 'width'."""
    ids = [e[0] for e in entries]
    dup_ids = [i for i, n in Counter(ids).items() if n > 1]
    if dup_ids:
        raise ValueError(f"{path}: barcode id(s) listed twice: {', '.join(dup_ids[:5])}")
    owner = {}
    for bc_id, seqs in entries:
        for s in seqs:
            if s in owner and owner[s] != bc_id:
                raise ValueError(f"{path}: sequence {s} belongs to both {owner[s]} and {bc_id}")
            owner[s] = bc_id
    longest = max(len(s) for _, seqs in entries for s in seqs)
    if longest > width:
        raise ValueError(f"{path}: segment '{name}' takes barcodes of at most {width} nt, found {longest} nt")
    # Pad shorter barcodes to the segment's width; pick a padding base that keeps every
    # canonical sequence unique
    for pad in PAD_BASES:
        canon = [seqs[0] + pad * (width - len(seqs[0])) for _, seqs in entries]
        if len(set(canon)) == len(canon):
            break
    else:
        raise ValueError(f"{path}: the barcodes of '{name}' are not unique once padded to {width} nt")
    return [
        {"id": bc_id, "sequences": seqs, "canonical": c}
        for (bc_id, seqs), c in zip(entries, canon)
    ]


def cmd_prepare(args):
    lay = LAYOUTS[args.layout]
    if len(args.whitelists) != len(lay["segments"]):
        sys.exit(
            f"Layout '{args.layout}' has {len(lay['segments'])} barcode segments "
            f"({' '.join(lay['segments'])}), but {len(args.whitelists)} barcode list(s) were given."
        )
    os.makedirs(args.outdir, exist_ok=True)
    out = {"layout": args.layout, "segments": []}
    for n, (name, path) in enumerate(zip(lay["segments"], args.whitelists), 1):
        try:
            entries = prepare_segment(name, read_barcode_list(path), lay["width"][name], path)
        except ValueError as exc:
            sys.exit(str(exc))
        out["segments"].append({"name": name, "source": os.path.basename(path), "entries": entries})
        with open(os.path.join(args.outdir, f"whitelist_{n}_{name}.txt"), "w") as fh:
            fh.writelines(e["canonical"] + "\n" for e in entries)
        print(f"{name}: {len(entries)} barcodes from {path}", file=sys.stderr)
    if lay.get("cb_code"):
        combos = 1
        for seg in out["segments"]:
            combos *= len(seg["entries"])
        if combos > 4 ** lay["cb_code"]:
            sys.exit(f"The barcode lists make {combos} cells, more than the {lay['cb_code']}-nt cell barcode holds ({4 ** lay['cb_code']}).")
    with open(os.path.join(args.outdir, "barcodes.json"), "w") as fh:
        json.dump(out, fh, indent=1)


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def build_match_table(seqs_by_index, max_mismatches):
    """
    {observed: (entry index, mismatches)} for every sequence within 'max_mismatches'
    (Hamming, over A/C/G/T/N) of exactly one entry at the smallest distance; observed
    sequences equally close to two entries map to (None, distance), i.e. ambiguous.
    """
    best = {}
    for idx, seqs in seqs_by_index:
        for seq in seqs:
            for k in range(max_mismatches + 1):
                for positions in itertools.combinations(range(len(seq)), k):
                    choices = [[b for b in "ACGTN" if b != seq[p]] for p in positions]
                    for subs in itertools.product(*choices):
                        obs = list(seq)
                        for p, b in zip(positions, subs):
                            obs[p] = b
                        obs = "".join(obs)
                        cur = best.get(obs)
                        if cur is None or k < cur[1]:
                            best[obs] = (idx, k)
                        elif k == cur[1] and cur[0] != idx:
                            best[obs] = (None, k)
    return best


class Segment:
    def __init__(self, seg, max_mismatches, revcomp_seqs=False):
        self.name = seg["name"]
        self.ids = [e["id"] for e in seg["entries"]]
        self.canonical = [e["canonical"] for e in seg["entries"]]
        seqs = [(i, [revcomp(s) if revcomp_seqs else s for s in e["sequences"]])
                for i, e in enumerate(seg["entries"])]
        self.table = build_match_table(seqs, max_mismatches)
        self.id_index = {bc_id: i for i, bc_id in enumerate(self.ids)}
        self.lengths = sorted({len(s) for e in seg["entries"] for s in e["sequences"]})
        self.stats = Counter()

    def match_within(self, seq):
        """
        A barcode at the end (or else the start) of a longer sequence, when its place in the
        tag is not given: the Ultima trimmer's segments may hold a constant linker next to the
        barcode (e.g. QuantumScale b2/b3: 4 nt, then 8 nt). Giving the place ('bead2=b2:4')
        is safer: a linker plus part of the barcode can come close to another barcode.
        """
        hits = set()
        for end in ("end", "start"):
            for length in self.lengths:
                if len(seq) > length:
                    hit = self.table.get(seq[-length:] if end == "end" else seq[:length])
                    if hit is not None:
                        hits.add(hit)
            if hits:
                break
        if not hits:
            self.stats["NoMatch"] += 1
            return None
        best = min(k for _, k in hits)
        idx = {i for i, k in hits if k == best}
        if len(idx) != 1 or None in idx:
            self.stats["Ambiguous"] += 1
            return None
        self.stats["Exact" if best == 0 else "Corrected"] += 1
        return idx.pop()

    def match_value(self, value):
        """
        A barcode from a SAM tag (Ultima): a barcode id of the list, a sequence, or a value
        holding one of these between separators (e.g. a read group 'X-QSR-1' or 'ID_ACGT').
        """
        idx = self.id_index.get(value)
        if idx is not None:
            self.stats["Exact"] += 1
            return idx
        seq = value.upper()
        if re.fullmatch(r"[ACGTN]+", seq):
            if len(seq) > self.lengths[-1]:
                return self.match_within(seq)
            return self.match(seq)
        parts = [p for p in re.split(r"[-_:.|/\s]+", value) if p]
        # ids may hold separators themselves ('QSR-1'): try every run of parts, longest first
        for n in range(len(parts), 0, -1):
            for i in range(len(parts) - n + 1):
                for sep in "-_":
                    idx = self.id_index.get(sep.join(parts[i:i + n]))
                    if idx is not None:
                        self.stats["Exact"] += 1
                        return idx
        for p in parts:
            if re.fullmatch(r"[ACGTNacgtn]+", p) and p.upper() in self.table:
                return self.match(p.upper())
        self.stats["NoMatch"] += 1
        return None

    def match(self, obs):
        """Entry index, or None; counts Exact/Corrected/Ambiguous/NoMatch."""
        hit = self.table.get(obs)
        if hit is None:
            self.stats["NoMatch"] += 1
            return None
        idx, k = hit
        if idx is None:
            self.stats["Ambiguous"] += 1
            return None
        self.stats["Exact" if k == 0 else "Corrected"] += 1
        return idx


# ---------------------------------------------------------------------------
# Wells and samples
# ---------------------------------------------------------------------------

WELL_RE = re.compile(r"^(\d{1,2})([A-Pa-p])$")


def _well_order_key(well_id):
    """Column-major plate order: 1A, 1B, ... 1H, 2A, ..."""
    m = WELL_RE.match(well_id)
    return (int(m.group(1)), m.group(2).upper())


def expand_wells(spec, ids):
    """
    Barcode ids named by a well specification: ids separated by ';' (or whitespace),
    each a single id or an inclusive range 'first-last'. Ranges of plate wells ('1A-2C')
    run column by column (1A..1H, 2A..2C); ranges of other ids follow the barcode list.
    An empty specification names every id.
    """
    spec = (spec or "").strip()
    if not spec:
        return list(ids)
    known = set(ids)
    out = []
    for part in re.split(r"[;\s]+", spec):
        if not part:
            continue
        if part in known:
            out.append(part)
            continue
        # Ids may hold dashes themselves (e.g. 'QSR-1'): split at the dash that leaves an id on both sides
        splits = [(part[:i], part[i + 1:]) for i, c in enumerate(part) if c == "-"]
        splits = [(a, b) for a, b in splits if a in known and b in known]
        if len(splits) != 1:
            raise ValueError(f"'{part}' is not a barcode id or a range of them")
        first, last = splits[0]
        if WELL_RE.match(first) and WELL_RE.match(last):
            lo, hi = _well_order_key(first), _well_order_key(last)
            if lo > hi:
                raise ValueError(f"range '{part}' runs backwards")
            out.extend(i for i in sorted((i for i in ids if WELL_RE.match(i)), key=_well_order_key)
                       if lo <= _well_order_key(i) <= hi)
        else:
            a, b = ids.index(first), ids.index(last)
            if a > b:
                raise ValueError(f"range '{part}' runs backwards")
            out.extend(ids[a:b + 1])
    seen = set()
    return [i for i in out if not (i in seen or seen.add(i))]


def read_samples(path, rt_ids, pcr_ids):
    """
    samples.tsv: 'sample<TAB>rt wells<TAB>pcr wells' per line. Returns the sample names
    and {(rt index, pcr index): sample}.
    """
    names, assign = [], {}
    rt_pos = {i: n for n, i in enumerate(rt_ids)}
    pcr_pos = {i: n for n, i in enumerate(pcr_ids)}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t") + ["", ""]
            sample, rt_spec, pcr_spec = fields[0].strip(), fields[1], fields[2]
            try:
                rts = expand_wells(rt_spec, rt_ids)
                pcrs = expand_wells(pcr_spec, pcr_ids)
            except ValueError as exc:
                sys.exit(f"Sample '{sample}': {exc}")
            names.append(sample)
            for r in rts:
                for p in pcrs:
                    key = (rt_pos[r], pcr_pos[p])
                    if key in assign and assign[key] != sample:
                        sys.exit(f"RT well {r} with PCR well {p} is given to both '{assign[key]}' and '{sample}'.")
                    assign[key] = sample
    if len(set(names)) != len(names):
        sys.exit(f"{path}: a sample is listed twice")
    return names, assign


# ---------------------------------------------------------------------------
# FASTQ I/O
# ---------------------------------------------------------------------------

def open_in(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    if gz and shutil.which("pigz"):
        proc = subprocess.Popen(["pigz", "-dc", path], stdout=subprocess.PIPE)
        return _ProcReader(proc)
    return gzip.open(path, "rt") if gz else open(path)


class _ProcReader:
    def __init__(self, proc):
        self.proc = proc
        self.fh = (line.decode() for line in proc.stdout)

    def __iter__(self):
        return self.fh

    def close(self):
        self.proc.stdout.close()
        self.proc.wait()


class GzWriter:
    """Gzipped output through pigz when it is there, else Python's gzip."""

    def __init__(self, path, threads=1):
        self.path = path
        if shutil.which("pigz"):
            self.out = open(path, "wb")
            self.proc = subprocess.Popen(["pigz", "-c", "-3", "-p", str(threads)],
                                         stdin=subprocess.PIPE, stdout=self.out)
            self.fh = self.proc.stdin
            self.encode = True
        else:
            self.proc = None
            self.fh = gzip.open(path, "wt", compresslevel=3)
            self.encode = False

    def write(self, text):
        self.fh.write(text.encode() if self.encode else text)

    def close(self):
        self.fh.close()
        if self.proc:
            if self.proc.wait() != 0:
                sys.exit(f"pigz failed writing {self.path}")
            self.out.close()


def fastq_records(path):
    fh = open_in(path)
    try:
        it = iter(fh)
        for header in it:
            seq = next(it).rstrip("\n")
            next(it)
            qual = next(it).rstrip("\n")
            yield header.rstrip("\n"), seq, qual
    finally:
        fh.close()


def chain_records(paths):
    for p in paths:
        yield from fastq_records(p)


# ---------------------------------------------------------------------------
# Trimming
# ---------------------------------------------------------------------------

def trim_cdna(seq, trim):
    """
    Cut the cDNA read at the first homopolymer run of a trimmed base at least as long as
    its threshold (e.g. a poly-A tail of 8 or more), then drop a run of 3 or more of a
    trimmed base left at the end of the read. Returns the kept length.
    """
    cut = len(seq)
    for base, n in trim.items():
        pos = seq.find(base * n)
        if pos != -1 and pos < cut:
            cut = pos
    changed = True
    while changed and cut:
        changed = False
        for base in trim:
            run = len(seq[:cut]) - len(seq[:cut].rstrip(base))
            if run >= 3:
                cut -= run
                changed = True
    return cut


# ---------------------------------------------------------------------------
# demux
# ---------------------------------------------------------------------------

READ_CATEGORIES = ["Pass", "LinkerError", "BarcodeError", "UmiFiltered", "TooShort", "Unassigned"]


def find_linker(seq, linker):
    """Start of the splint in the barcode read, or None."""
    ref = linker["sequence"]
    for k in range(linker["max_mismatches"] + 1):
        for s in linker["starts"]:
            obs = seq[s:s + len(ref)]
            if len(obs) == len(ref) and sum(a != b for a, b in zip(obs, ref)) == k:
                return s
    return None


def umi_ok(umi):
    return "N" not in umi and len(set(umi)) > 1


def fastq_units(reads, lay):
    """
    Per read of a FASTQ chunk: (name, cDNA seq, cDNA qual, error category or None,
    {segment: observed barcode or None}, UMI).
    """
    linker = lay["linker"]
    umi_src, umi_off, umi_len = lay["umi"]
    extract = lay["extract"]
    iters = [chain_records(reads["cdna"]), chain_records(reads["bc"])]
    iters += [chain_records(reads[r]) for r in ("i1", "i2") if reads[r]]
    has_i1, has_i2 = bool(reads["i1"]), bool(reads["i2"])
    for recs in zip(*iters):
        header, cseq, cqual = recs[0]
        name = header.split(None, 1)[0]
        bseq = recs[1][1]
        seqs = {"bc": bseq}
        if has_i1:
            seqs["i1"] = recs[2][1]
        if has_i2:
            seqs["i2"] = recs[-1][1]

        anchor = None
        if linker:
            anchor = find_linker(bseq, linker)
            if anchor is None:
                yield name, cseq, cqual, "LinkerError", None, None
                continue

        obs = {}
        for seg, (src, off, length) in extract.items():
            if src == "anchor":
                o = bseq[anchor + off:anchor + off + length]
            elif length is None:      # up to the linker
                o = seqs[src][off:anchor]
            else:
                o = seqs[src][off:off + length]
            obs[seg] = o if length is None or len(o) == length else None

        if umi_src == "anchor":
            umi = bseq[anchor + umi_off:anchor + umi_off + umi_len]
        else:
            umi = seqs[umi_src][umi_off:umi_off + umi_len]
        yield name, cseq, cqual, None, obs, umi

    # A read stream that ends early would silently drop the rest of the longer ones
    for it in iters:
        if next(it, None) is not None:
            sys.exit("The read files of this chunk hold different numbers of reads.")


def parse_ultima_umi(spec):
    """'tag:XX' (a SAM tag) or 'seq:START:LENGTH' (cut from the read, 0-based)."""
    m = re.fullmatch(r"tag:([A-Za-z][A-Za-z0-9])", spec or "")
    if m:
        return ("tag", m.group(1))
    m = re.fullmatch(r"seq:(\d+):(\d+)", spec or "")
    if m:
        return ("seq", int(m.group(1)), int(m.group(2)))
    sys.exit(f"Ultima UMI location '{spec}': give 'tag:XX' (the SAM tag holding it) or 'seq:START:LENGTH' "
             "(its place in the read, 0-based). Look at a read with 'samtools view FILE.cram | head -1'.")


def parse_ultima_tags(spec, segments):
    """
    'segment=TAG[:START],...' for every barcode segment: {segment: (TAG, START or None)}.
    START is where the barcode begins in the tag's value (0-based), for tags that hold a
    constant linker before it.
    """
    tags = {}
    for part in re.split(r"[,\s]+", (spec or "").strip()):
        if not part:
            continue
        m = re.fullmatch(r"(\w+)=([A-Za-z][A-Za-z0-9])(?::(\d+))?", part)
        if not m:
            sys.exit(f"Ultima tags: '{part}' is not 'segment=TAG' or 'segment=TAG:START'")
        tags[m.group(1)] = (m.group(2), int(m.group(3)) if m.group(3) else None)
    missing = [s for s in segments if s not in tags]
    if missing:
        sys.exit(f"Ultima tags: no tag given for {', '.join(missing)} (got '{spec}')")
    return tags


def sam_records(paths, reference=None, threads=1):
    """(name, flag, seq, qual, {tag: value}) of every primary record of SAM/BAM/CRAM files."""
    for path in paths:
        cmd = ["samtools", "view", "-@", str(threads)]
        if reference:
            cmd += ["-T", reference]
        proc = subprocess.Popen(cmd + [path], stdout=subprocess.PIPE, text=True)
        finished = False
        try:
            for line in proc.stdout:
                f = line.rstrip("\n").split("\t")
                flag = int(f[1])
                if flag & 0x900:      # secondary, supplementary
                    continue
                seq, qual = f[9], f[10]
                if flag & 0x10:
                    seq, qual = revcomp(seq), qual[::-1]
                tags = {t[:2]: t[5:] for t in f[11:]}
                yield f[0], flag, seq, qual, tags
            finished = True
        finally:
            if not finished:      # the reader stopped early (e.g. the orientation check)
                proc.kill()
            proc.stdout.close()
            if proc.wait() != 0 and finished:
                sys.exit(f"samtools view failed on {path}")


def ultima_units(paths, seg_tags, umi_spec, reference=None, threads=1):
    """Per record of Ultima CRAMs, whose barcodes the Ultima trimmer moved into tags."""
    for name, _flag, seq, qual, tags in sam_records(paths, reference, threads):
        obs = {}
        for seg, (tag, start) in seg_tags.items():
            value = tags.get(tag)
            obs[seg] = value[start:] if value is not None and start else value
        if umi_spec[0] == "tag":
            umi = tags.get(umi_spec[1], "")
        else:
            start, length = umi_spec[1], umi_spec[2]
            umi = seq[start:start + length]
            seq, qual = seq[:start] + seq[start + length:], qual[:start] + qual[start + length:]
        yield "@" + name, seq, qual, None, obs, umi


def detect_pcr_revcomp(observed, pcr_seg_entries, n_reads=200000):
    """True when more of the first observed PCR barcodes match the list reverse-complemented."""
    fwd = {s for e in pcr_seg_entries for s in e["sequences"]}
    rev = {revcomp(s) for s in fwd}
    n_fwd = n_rev = 0
    for n, obs in enumerate(observed):
        if n >= n_reads:
            break
        n_fwd += obs in fwd
        n_rev += obs in rev
    return n_rev > n_fwd, n_fwd, n_rev


def cmd_demux(args):
    with open(args.barcodes) as fh:
        bcs = json.load(fh)
    layout = bcs["layout"]
    lay = LAYOUTS[layout]
    segs_json = {s["name"]: s for s in bcs["segments"]}
    seg_order = lay["segments"]

    ultima = bool(args.cram)
    if ultima:
        if args.cdna or args.bc or args.i1 or args.i2:
            sys.exit("Give either --cram (Ultima) or the FASTQs, not both.")
        seg_tags = parse_ultima_tags(args.ultima_tags, seg_order)
        umi_spec = parse_ultima_umi(args.ultima_umi)
    else:
        reads = {"cdna": args.cdna or [], "bc": args.bc or [], "i1": args.i1 or [], "i2": args.i2 or []}
        for r in ["cdna", "bc"] + lay["index_reads"]:
            if not reads[r]:
                sys.exit(f"Layout '{layout}' needs the {r.upper()} read, which was not given.")
        n_files = {r: len(reads[r]) for r in ["cdna", "bc"] + lay["index_reads"]}
        if len(set(n_files.values())) != 1:
            sys.exit(f"Every read needs the same number of files: {n_files}")

    # The i5 index is read reverse-complemented on some sequencers
    i5_rc_detect = None
    if args.i5_revcomp == "auto":
        if ultima:
            observed = (tags.get(seg_tags["pcr"][0], "")[seg_tags["pcr"][1] or 0:]
                        for _, _, _, _, tags in sam_records(args.cram, args.reference, args.threads))
        else:
            _, off, length = lay["extract"]["pcr"]
            observed = (seq[off:off + length] for _, seq, _ in chain_records(reads["i2"]))
        rc, n_fwd, n_rev = detect_pcr_revcomp(observed, segs_json["pcr"]["entries"])
        i5_rc_detect = {"forward_exact": n_fwd, "revcomp_exact": n_rev}
    else:
        rc = args.i5_revcomp == "true"
    print(f"PCR index read {'reverse-complemented' if rc else 'as listed'}", file=sys.stderr)

    segments = {name: Segment(segs_json[name], args.max_mismatches, revcomp_seqs=(name == "pcr" and rc))
                for name in seg_order}
    names, assign = read_samples(args.samples, segments["rt"].ids, segments["pcr"].ids)
    rt_idx_of = seg_order.index("rt")
    pcr_idx_of = seg_order.index("pcr")

    os.makedirs(args.outdir, exist_ok=True)
    writers = {}
    for s in names + (["unassigned"] if args.write_unassigned else []):
        writers[s] = (GzWriter(os.path.join(args.outdir, f"{s}.{args.chunk}_cDNA.fastq.gz")),
                      GzWriter(os.path.join(args.outdir, f"{s}.{args.chunk}_BC_UMI.fastq.gz")))

    categories = Counter()
    per_sample = Counter()
    per_rt = Counter()
    umi_len = lay["umi"][2]
    trim = lay["trim"]
    min_len = args.min_cdna_length
    syn_qual = "F" * synthetic_length(layout)
    cb_code = lay.get("cb_code")
    sizes = [len(segments[n].ids) for n in seg_order]

    if ultima:
        units = ultima_units(args.cram, seg_tags, umi_spec, args.reference, args.threads)
    else:
        units = fastq_units(reads, lay)
    total = 0
    for name, cseq, cqual, error, obs, umi in units:
        total += 1
        if error:
            categories[error] += 1
            continue

        hits = []
        for seg in seg_order:
            o = obs[seg]
            if o is None:
                segments[seg].stats["Error"] += 1
                hits.append(None)
            else:
                hits.append(segments[seg].match_value(o) if ultima else segments[seg].match(o))
        if None in hits:
            categories["BarcodeError"] += 1
            continue

        if len(umi) != umi_len or not umi_ok(umi):
            categories["UmiFiltered"] += 1
            continue

        keep = trim_cdna(cseq, trim)
        if keep < min_len:
            categories["TooShort"] += 1
            continue

        rt_hit = hits[rt_idx_of]
        per_rt[segments["rt"].ids[rt_hit]] += 1
        sample = assign.get((rt_hit, hits[pcr_idx_of]))
        if sample is None:
            categories["Unassigned"] += 1
            if not args.write_unassigned:
                continue
            sample = "unassigned"
        else:
            categories["Pass"] += 1
            per_sample[sample] += 1

        if cb_code:
            cb = encode_cell(hits, sizes, cb_code)
        else:
            cb = "".join(segments[n].canonical[h] for n, h in zip(seg_order, hits))
        w_cdna, w_bc = writers[sample]
        w_cdna.write(f"{name}\n{cseq[:keep]}\n+\n{cqual[:keep]}\n")
        w_bc.write(f"{name}\n{cb}{umi}\n+\n{syn_qual}\n")

    for w_cdna, w_bc in writers.values():
        w_cdna.close()
        w_bc.close()

    metrics = {
        "chunk": args.chunk,
        "layout": layout,
        "input": "ultima_cram" if ultima else "fastq",
        "i5_revcomp": rc,
        "i5_orientation_check": i5_rc_detect,
        "reads": total,
        "read_categories": {c: categories[c] for c in READ_CATEGORIES},
        "barcodes": {n: dict(segments[n].stats) for n in seg_order},
        "samples": {s: per_sample[s] for s in names},
        "rt_wells": dict(sorted(per_rt.items())),
    }
    with open(os.path.join(args.outdir, f"{args.chunk}.metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=1)
    print(json.dumps(metrics["read_categories"]), file=sys.stderr)


# ---------------------------------------------------------------------------
# gather
# ---------------------------------------------------------------------------

def cmd_gather(args):
    """One sample's metrics: the sum of every chunk of the libraries it was read from."""
    reads = 0
    cats, rt = Counter(), Counter()
    barcodes = {}
    chunks, sample_reads = [], 0
    for path in sorted(args.metrics):
        with open(path) as fh:
            m = json.load(fh)
        chunks.append({"chunk": m["chunk"], "i5_revcomp": m["i5_revcomp"], "reads": m["reads"]})
        reads += m["reads"]
        cats.update(m["read_categories"])
        rt.update(m["rt_wells"])
        for seg, st in m["barcodes"].items():
            barcodes.setdefault(seg, Counter()).update(st)
        sample_reads += m["samples"].get(args.sample, 0)
    out = {
        "sample": args.sample,
        "sample_reads": sample_reads,
        "library_reads": reads,
        "library_read_categories": {c: cats[c] for c in READ_CATEGORIES},
        "library_barcodes": {s: dict(c) for s, c in barcodes.items()},
        "library_rt_wells": dict(sorted(rt.items())),
        "chunks": chunks,
    }
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)


# ---------------------------------------------------------------------------
# decode
# ---------------------------------------------------------------------------

def cmd_decode(args):
    """Cell barcodes (one per line, e.g. a barcodes.tsv) back to their barcode ids per segment."""
    with open(args.barcodes) as fh:
        bcs = json.load(fh)
    lay = LAYOUTS[bcs["layout"]]
    segs = bcs["segments"]
    sizes = [len(s["entries"]) for s in segs]
    by_canon = [{e["canonical"]: e["id"] for e in s["entries"]} for s in segs]
    out = sys.stdout
    out.write("barcode\t" + "\t".join(s["name"] for s in segs) + "\n")
    with (gzip.open if args.cells.endswith(".gz") else open)(args.cells, "rt") as fh:
        for line in fh:
            cb = line.strip().split("\t")[0]
            if not cb:
                continue
            ids = []
            if lay.get("cb_code"):
                try:
                    ids = [s["entries"][i]["id"] for s, i in zip(segs, decode_cell(cb, sizes))]
                except ValueError:
                    ids = []
            else:
                pos = 0
                for name, lookup in zip(lay["segments"], by_canon):
                    w = lay["width"][name]
                    ids.append(lookup.get(cb[pos:pos + w], "?"))
                    pos += w
            out.write(cb + "\t" + "\t".join(ids or ["?"] * len(segs)) + "\n")


# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prepare", help="check the barcode lists and write the mapping whitelists")
    p.add_argument("--layout", required=True, choices=sorted(LAYOUTS))
    p.add_argument("--whitelists", nargs="+", required=True, help="one barcode list per segment, in segment order")
    p.add_argument("--outdir", default=".")
    p.set_defaults(func=cmd_prepare)

    p = sub.add_parser("demux", help="correct barcodes and split one chunk by sample")
    p.add_argument("--barcodes", required=True, help="barcodes.json from 'prepare'")
    p.add_argument("--samples", required=True, help="TSV: sample, RT wells, PCR wells (empty = all)")
    p.add_argument("--chunk", required=True, help="name of this chunk, used in the output names")
    p.add_argument("--cdna", nargs="+", help="cDNA read file(s)")
    p.add_argument("--bc", nargs="+", help="barcode/UMI read file(s)")
    p.add_argument("--i1", nargs="*", help="index 1 read file(s)")
    p.add_argument("--i2", nargs="*", help="index 2 read file(s)")
    p.add_argument("--cram", nargs="+", help="Ultima CRAM/BAM file(s) instead of FASTQs: barcodes in tags, cDNA in the read")
    p.add_argument("--ultima-tags", default="pcr=RG,bead1=b1,bead2=b2:4,bead3=b3:4,rt=b4",
                   help="with --cram: the tag of every barcode segment, as 'segment=TAG[:START],...' (START: where the barcode begins in the value)")
    p.add_argument("--ultima-umi", help="with --cram: where the UMI is, 'tag:XX' or 'seq:START:LENGTH' (0-based, cut from the read)")
    p.add_argument("--reference", help="with --cram: FASTA for reference-compressed CRAMs")
    p.add_argument("--threads", type=int, default=1, help="with --cram: samtools decompression threads")
    p.add_argument("--max-mismatches", type=int, default=1, choices=[0, 1, 2])
    p.add_argument("--i5-revcomp", default="auto", choices=["auto", "true", "false"])
    p.add_argument("--min-cdna-length", type=int, default=16)
    p.add_argument("--write-unassigned", action="store_true")
    p.add_argument("--outdir", default=".")
    p.set_defaults(func=cmd_demux)

    p = sub.add_parser("gather", help="sum the chunk metrics for one sample")
    p.add_argument("--sample", required=True)
    p.add_argument("--metrics", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_gather)

    p = sub.add_parser("decode", help="cell barcodes back to their barcode ids per segment")
    p.add_argument("--barcodes", required=True, help="barcodes.json from 'prepare'")
    p.add_argument("--cells", required=True, help="one cell barcode per line (first column), plain or gzipped")
    p.set_defaults(func=cmd_decode)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
