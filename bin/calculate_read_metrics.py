#!/usr/bin/env python3
"""
Read-level mtDNA, rRNA and antisense metrics of a STARsolo run (CALC_READ_METRICS).

Writes two key,value CSV files and one per-barcode table:

``--out-mt-rrna``
    mtDNA and rRNA reads, per library and per called cell, from ONE streaming
    pass over the BAM. Every category -- unique, multimapped, primary alignment,
    all alignments, called cells, mitochondrial reads by origin -- is tallied in
    that pass, and rRNA / mitochondrial-gene overlaps are tested in-process on
    the aligned blocks of each record, the rule featureCounts uses. Nothing is
    written to disk but the metrics.

``--out-barcode-reads``
    The same pass's mapped, mtDNA and rRNA reads (primary alignment) for every
    barcode, as a gzipped TSV. PERCELL_METRICS reads it instead of scanning the
    BAM again, so the per-cell and the called-cell values share one count.

``--out-antisense``
    The antisense share of gene-overlapping reads, from STARsolo's
    ``CellReads.stats`` rather than from the BAM. STARsolo flags every read with a
    unique genomic alignment as exonic, intronic, exonicAS or intronicAS against
    the same GeneFull_Ex50pAS model and strand it counts the matrix with, so the
    numbers describe exactly the reads behind the matrix, introns included.
    Written only for a stranded run (Forward or Reverse).

Metric names carry commas ("of mapped reads, primary alignment"), so they are
quoted; readers must parse the files as CSV.
"""

import argparse
import csv
import gzip
import os
import re
import sys
from bisect import bisect_right
from collections import Counter, defaultdict, namedtuple

import pysam


# ------------------------------------------------------------------
# GTF conventions
# ------------------------------------------------------------------
# Ensembl spells the biotype gene_biotype on gene rows; GENCODE and NCBI-derived
# GTFs use gene_type, and some annotations carry it only on transcript or exon
# rows. Listed in order of preference, which breaks ties.
FEATURE_TYPES = ("gene", "transcript", "exon")
BIOTYPE_ATTRS = ("gene_biotype", "gene_type", "transcript_biotype", "transcript_type")

# Every spelling at once: a concatenation of ref_gtf and ref_gtf_addfeature may mix them
_ANY_BIOTYPE_RE = re.compile(r'(?:gene_biotype|gene_type|transcript_biotype|transcript_type) "([^"]*)"')
_ATTR_VALUE_RE = {a: re.compile(a + r' "([^"]*)"', re.IGNORECASE) for a in BIOTYPE_ATTRS}

# BAM flags: a "read" is counted once, at its primary alignment. STAR writes one
# record per alignment, so a read with NH:i:6 has six records; secondary (0x100)
# and supplementary (0x800) records are skipped wherever reads are counted. Only
# the rows labelled "alignments" count every record, on purpose.
FLAG_UNMAPPED = 0x4
FLAG_NOT_PRIMARY = 0x100 | 0x800

# Read categories of the rRNA and mtDNA rows: (Counts key, unit, pool, qualifier).
# mapped/unique/mmpa count reads at their primary alignment; mmaa every alignment
# of every multimapper. The called-cell rows use the first two.
CATEGORIES = (
    ("mapped", "reads", "mapped reads", "primary alignment"),
    ("unique", "reads", "uniquely mapped reads", None),
    ("mmpa", "reads", "multimapped reads", "primary alignment"),
    ("mmaa", "alignments", "multimapped read alignments", "all alignments"),
)
CELL_CATEGORIES = CATEGORIES[:2]


def warn(*lines):
    for line in lines:
        print(f"WARNING: {line}", file=sys.stderr)


def usable_file(path, what):
    """True when *path* is a non-empty file; warns about it otherwise."""
    if os.path.isfile(path) and os.path.getsize(path) > 0:
        return True
    warn(f"{what} is empty or unreadable: {path}")
    return False


# ------------------------------------------------------------------
# Intervals
# ------------------------------------------------------------------

def merge_intervals(intervals):
    """Sorted, non-overlapping ``(starts, ends)`` lists (0-based, half-open) for bisection."""
    starts, ends = [], []
    for start, end in sorted(intervals):
        if ends and start <= ends[-1]:
            ends[-1] = max(ends[-1], end)
        else:
            starts.append(start)
            ends.append(end)
    return starts, ends


def overlaps(blocks, region):
    """True when any aligned block overlaps an interval of *region* by >= 1 bp.

    Aligned blocks rather than the reference span, as featureCounts does: a read
    spliced across an intron does not count for a feature inside that intron.
    """
    if region is None:
        return False
    starts, ends = region
    n = len(starts)
    for bstart, bend in blocks:
        i = bisect_right(starts, bstart) - 1
        if i >= 0 and ends[i] > bstart:
            return True
        if i + 1 < n and starts[i + 1] < bend:
            return True
    return False


def by_tid(intervals_by_contig, tid_of):
    """``{contig: [intervals]}`` to ``{tid: (starts, ends)}`` for the contigs in the BAM."""
    return {tid_of[c]: merge_intervals(ivs) for c, ivs in intervals_by_contig.items()
            if ivs and c in tid_of}


# ------------------------------------------------------------------
# Annotation
# ------------------------------------------------------------------

def collect_rrna_biotypes(line, rrna_values):
    """Add the biotype values on *line* containing "rRNA" to ``rrna_values[attr]``.

    Collected from every line, under each spelling, case-insensitively; only the
    winning spelling's are reported.
    """
    if "rrna" not in line.lower():
        return
    for attr, regex in _ATTR_VALUE_RE.items():
        for value in regex.findall(line):
            if "rrna" in value.lower():
                rrna_values[attr].add(value)


def best_biotype_pair(pair_rows):
    """``(feature_type, biotype_attr, rows)`` of the pair annotating the most rows.

    ``(None, None, 0)`` when no row carries a biotype; ties go to the pair listed
    first in FEATURE_TYPES and BIOTYPE_ATTRS.
    """
    best_type, best_attr, best = None, None, 0
    for ftype in FEATURE_TYPES:
        for attr in BIOTYPE_ATTRS:
            if pair_rows[(ftype, attr)] > best:
                best, best_type, best_attr = pair_rows[(ftype, attr)], ftype, attr
    return best_type, best_attr, best


def parse_gtf(gtf_path, mt_contigs):
    """Everything the metrics need from the main annotation, in one pass.

    Returns a dict with:
      feature_type, biotype_attr  the pair annotating the most rows (None, None if none)
      biotype_rows, feature_rows  how many rows of that type carry the attribute, of all
      rrna_rows                   [(contig, start0, end)] rRNA features of feature_type,
                                  under any biotype spelling
      rrna_biotypes               distinct biotype values containing "rRNA", "|"-joined
      mt_exons                    {"+"|"-": {contig: [intervals]}} exons on mt_contigs
      mt_exon_rows                number of those exon rows

    The winning pair is the one annotating the MOST rows, rather than the first
    one occurring anywhere: in a concatenation of ref_gtf and ref_gtf_addfeature,
    a handful of hand-written spike-in lines using gene_biotype must not outvote an
    entire GENCODE annotation using gene_type, or every rRNA gene of the base
    annotation would be skipped.
    """
    mt_set = set(mt_contigs)
    rows = defaultdict(int)
    pair_rows = defaultdict(int)
    rrna_by_type = {ft: [] for ft in FEATURE_TYPES}
    rrna_values = {a: set() for a in BIOTYPE_ATTRS}
    mt_exons = {"+": defaultdict(list), "-": defaultdict(list)}
    mt_exon_rows = 0

    with open(gtf_path) as gtf:
        for line in gtf:
            collect_rrna_biotypes(line, rrna_values)
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            contig, ftype, attributes = cols[0], cols[2], cols[8]

            if ftype in rrna_by_type:
                rows[ftype] += 1
                for attr in BIOTYPE_ATTRS:
                    if attr + ' "' in attributes:
                        pair_rows[(ftype, attr)] += 1
                if any("rrna" in v.lower() for v in _ANY_BIOTYPE_RE.findall(attributes)):
                    rrna_by_type[ftype].append((contig, int(cols[3]) - 1, int(cols[4])))

            if ftype == "exon" and contig in mt_set:
                strand = "-" if cols[6] == "-" else "+"
                mt_exons[strand][contig].append((int(cols[3]) - 1, int(cols[4])))
                mt_exon_rows += 1

    best_type, best_attr, best = best_biotype_pair(pair_rows)
    return {
        "feature_type": best_type,
        "biotype_attr": best_attr,
        "biotype_rows": best,
        "feature_rows": rows[best_type] if best_type else 0,
        "rrna_rows": rrna_by_type[best_type] if best_type else [],
        "rrna_biotypes": "|".join(sorted(rrna_values[best_attr])) if best_attr else "",
        "mt_exons": mt_exons,
        "mt_exon_rows": mt_exon_rows,
    }


def added_contigs(rrna_gtf_path, contig_lengths):
    """Contigs of the added rRNA reference, split into ``(present, missing)`` in the BAM.

    Those files are the added rRNA reference (ref_fasta_addfeature carries the
    sequences), so every read aligning to one of their contigs is an rRNA read,
    whether or not the added GTF spells out a biotype -- it commonly does not.
    """
    contigs = set()
    with open(rrna_gtf_path) as gtf:
        for line in gtf:
            cols = line.split("\t")
            if not line.startswith("#") and len(cols) >= 9:
                contigs.add(cols[0])
    ordered = sorted(contigs)
    return ([c for c in ordered if c in contig_lengths],
            [c for c in ordered if c not in contig_lengths])


def load_barcodes(path):
    """The called-cell barcodes: first column of a plain or gzipped barcodes.tsv."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return {line.split("\t")[0].strip() for line in fh if line.strip()}


def load_cells(path):
    """The called-cell barcodes, or None without a usable barcode list."""
    if path and usable_file(path, "called-cell barcode list"):
        return load_barcodes(path)
    return None


# ------------------------------------------------------------------
# Regions
# ------------------------------------------------------------------

# What the BAM pass counts as rRNA and mtDNA, keyed by BAM tid:
#   rrna_by_tid, mt_plus, mt_minus   {tid: (starts, ends)} for overlaps()
#   n_rrna_regions, added_present    what the rRNA regions were built from
#   mt_tids                          the mitochondrial contigs present in the BAM
#   has_mt_genes                     exon rows on those contigs, to split mtDNA by origin
#   sense_on_read_strand             True Forward, False Reverse, None unstranded
Regions = namedtuple("Regions", (
    "rrna_by_tid", "n_rrna_regions", "added_present",
    "mt_tids", "mt_plus", "mt_minus", "has_mt_genes", "sense_on_read_strand",
))


def present_mt_contigs(mt_contigs, tid_of):
    """The mitochondrial contigs found in the BAM header.

    A contig name that is not in the BAM header can only ever yield 0 reads,
    indistinguishable from a mitochondria-free library, so name the absent ones.
    """
    present = [c for c in mt_contigs if c in tid_of]
    if not mt_contigs:
        warn("no mitochondrial contig given; all mtDNA metrics will be 0")
    elif len(present) < len(mt_contigs):
        missing = " ".join(c for c in mt_contigs if c not in tid_of)
        warn(f"mitochondrial contig(s) absent from the BAM header: {missing}",
             "check mt_contig against the reference used for mapping")
    return present


def warn_biotypes(ann, gtf):
    """Warn when the rRNA biotypes of the annotation are partial or missing."""
    if ann["biotype_attr"] and ann["biotype_rows"] < ann["feature_rows"]:
        warn(f"{ann['biotype_attr']} annotates only {ann['biotype_rows']} of "
             f"{ann['feature_rows']} {ann['feature_type']} rows in {gtf}",
             "rows spelling the biotype differently are not counted -- typical when",
             "ref_gtf and ref_gtf_addfeature follow different GTF conventions")
    if not ann["rrna_biotypes"]:
        warn(f"no rRNA biotype found in {gtf}",
             f"biotype attribute detected: '{ann['biotype_attr'] or 'none'}', "
             f"feature type: '{ann['feature_type'] or 'none'}'")


def rrna_intervals(rrna_rows, rrna_gtf, lengths):
    """``({contig: [intervals]}, added_present)``: the rRNA regions to count.

    The rRNA features of the main annotation plus the whole of every contig of the
    added rRNA reference. Overlapping regions are merged later, so a read covered
    by both still counts once.
    """
    rrna = defaultdict(list)
    for contig, start, end in rrna_rows:
        rrna[contig].append((start, end))
    added_present = []
    if rrna_gtf and usable_file(rrna_gtf, "added rRNA annotation"):
        added_present, added_missing = added_contigs(rrna_gtf, lengths)
        for contig in added_present:
            rrna[contig].append((0, lengths[contig]))
        if added_missing:
            warn("contig(s) of the added rRNA reference are absent from the BAM header: "
                 + " ".join(added_missing),
                 "the reads were mapped against a reference built without ref_fasta_addfeature")
    return rrna, added_present


def mt_strand_sense(strand, has_mt_genes, gtf):
    """sense_on_read_strand for *strand*; warns when mtDNA reads cannot be split by origin."""
    if not has_mt_genes:
        warn(f"no exon rows on the mitochondrial contig(s) in {gtf};",
             "mtDNA reads are not split by origin")
    sense_on_read_strand = {"Forward": True, "Reverse": False}.get(strand)
    if has_mt_genes and sense_on_read_strand is None:
        warn(f"strand '{strand or 'unset'}' is not Forward or Reverse; mtDNA reads are",
             "not split into sense and antisense to mitochondrial genes")
    return sense_on_read_strand


def load_regions(args, bam):
    """``(annotation, Regions)`` for *bam*, from the GTFs and contigs in *args*."""
    tid_of = {name: i for i, name in enumerate(bam.references)}
    lengths = dict(zip(bam.references, bam.lengths))
    mt_present = present_mt_contigs(args.mt_contig, tid_of)

    ann = parse_gtf(args.gtf, mt_present)
    warn_biotypes(ann, args.gtf)
    rrna, added_present = rrna_intervals(ann["rrna_rows"], args.rrna_gtf, lengths)
    n_rrna_regions = len(ann["rrna_rows"]) + len(added_present)
    if n_rrna_regions == 0:
        warn(f"no rRNA regions to count: no rRNA biotype in {args.gtf} and no usable",
             "ref_gtf_addfeature contigs; all rRNA metrics will be reported as N/A")

    has_mt_genes = ann["mt_exon_rows"] > 0
    regions = Regions(
        rrna_by_tid=by_tid(rrna, tid_of),
        n_rrna_regions=n_rrna_regions,
        added_present=added_present,
        mt_tids={tid_of[m] for m in mt_present},
        mt_plus=by_tid(ann["mt_exons"]["+"], tid_of),
        mt_minus=by_tid(ann["mt_exons"]["-"], tid_of),
        has_mt_genes=has_mt_genes,
        sense_on_read_strand=mt_strand_sense(args.strand, has_mt_genes, args.gtf),
    )
    return ann, regions


# ------------------------------------------------------------------
# BAM pass
# ------------------------------------------------------------------

class Tally:
    """Reads (or alignments) of one category, and how many of them are rRNA and mtDNA."""

    __slots__ = ("total", "rrna", "mt")

    def __init__(self):
        self.total = self.rrna = self.mt = 0

    def add(self, is_rrna, is_mt):
        self.total += 1
        self.rrna += is_rrna
        self.mt += is_mt


class Counts:
    """Every tally of the BAM pass.

    library      {CATEGORIES key: Tally} over every record of the BAM
    in_cells     {CELL_CATEGORIES key: Tally} over the called-cell barcodes
    per_barcode  {CB: Tally} of primary alignments, for PERCELL_METRICS
    mt_*         mtDNA reads (primary alignment) by origin
    """

    def __init__(self):
        self.unmapped = 0
        self.library = {key: Tally() for key, *_ in CATEGORIES}
        self.in_cells = {key: Tally() for key, *_ in CELL_CATEGORIES}
        self.per_barcode = {}
        self.cells_found = set()
        self.mt_sense = self.mt_antisense = self.mt_in_genes = 0


def tag_or(read, name, default):
    """The value of tag *name* on *read*, or *default* when the record has none."""
    try:
        return read.get_tag(name)
    except KeyError:
        return default


def scan_bam(bam, regions, cells):
    """Tally every metric in one pass over *bam*, unmapped records included."""
    c = Counts()
    for read in bam.fetch(until_eof=True):
        count_record(c, read, regions, cells)
    return c


def count_record(c, read, regions, cells):
    """Tally one BAM record into *c*."""
    flag = read.flag
    if flag & FLAG_UNMAPPED:
        c.unmapped += 1
        return

    # A record without NH is neither unique nor multimapped, as before
    nh = tag_or(read, "NH", 0)
    primary = not flag & FLAG_NOT_PRIMARY
    multi = nh > 1
    if not primary and not multi:
        return

    tid = read.reference_id
    is_mt = tid in regions.mt_tids
    blocks = None
    is_rrna = False
    rrna_region = regions.rrna_by_tid.get(tid)
    if rrna_region is not None:
        blocks = read.get_blocks()
        is_rrna = overlaps(blocks, rrna_region)

    library = c.library
    if multi:
        # Every alignment of every multimapper, primary or not
        library["mmaa"].add(is_rrna, is_mt)
    if not primary:
        return

    unique = nh == 1
    library["mapped"].add(is_rrna, is_mt)
    if unique:
        library["unique"].add(is_rrna, is_mt)
    elif multi:
        library["mmpa"].add(is_rrna, is_mt)

    if is_mt and regions.has_mt_genes:
        count_mt_origin(c, read, tid, blocks, regions)
    cb = tag_or(read, "CB", None)
    if cb is not None:
        count_barcode(c, cb, is_rrna, is_mt, unique, cells)


def count_mt_origin(c, read, tid, blocks, regions):
    """Tally an mtDNA read as sense, antisense or outside the mitochondrial genes."""
    if blocks is None:
        blocks = read.get_blocks()
    on_plus = overlaps(blocks, regions.mt_plus.get(tid))
    on_minus = overlaps(blocks, regions.mt_minus.get(tid))
    if on_plus or on_minus:
        c.mt_in_genes += 1
    if regions.sense_on_read_strand is None:
        return
    # Genes on the read's own strand are sense in a Forward library, antisense in a Reverse one
    same, opposite = (on_minus, on_plus) if read.is_reverse else (on_plus, on_minus)
    if not regions.sense_on_read_strand:
        same, opposite = opposite, same
    c.mt_sense += same
    c.mt_antisense += opposite


def count_barcode(c, cb, is_rrna, is_mt, unique, cells):
    """Tally a primary alignment for its barcode, matrix or not, and for the called cells."""
    tally = c.per_barcode.get(cb)
    if tally is None:
        tally = c.per_barcode[cb] = Tally()
    tally.add(is_rrna, is_mt)

    if cells is not None and cb in cells:
        c.cells_found.add(cb)
        c.in_cells["mapped"].add(is_rrna, is_mt)
        if unique:
            c.in_cells["unique"].add(is_rrna, is_mt)


# ------------------------------------------------------------------
# Output
# ------------------------------------------------------------------

NA = "N/A"


def frac(num, den):
    """A 0-1 fraction, or N/A for a missing numerator or an empty denominator."""
    if num == NA or den == NA or not den:
        return NA
    return f"{num / den:.4f}"


class MetricsWriter:
    """``"name",value`` rows; free-text values are quoted as well."""

    def __init__(self, path):
        self.fh = open(path, "w", newline="")
        self.fh.write("Metric,Count\n")

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def row(self, name, value):
        self.fh.write(f'"{name}",{value}\n')

    def text(self, name, value):
        self.fh.write(f'"{name}","{value}"\n')

    def close(self):
        self.fh.close()


def write_barcode_reads(path, per_barcode):
    """Mapped, mtDNA and rRNA reads per barcode, as a gzipped TSV for PERCELL_METRICS.

    Every CB seen on a primary alignment, whether or not it is in the matrix; the
    reader keeps the barcodes it needs. Most reads first, so the file opens on cells.
    """
    with gzip.open(path, "wt", newline="") as fh:
        fh.write("CB\tMappedReads\tMTReads\trRNAReads\n")
        for cb, t in sorted(per_barcode.items(), key=lambda kv: -kv[1].total):
            fh.write(f"{cb}\t{t.total}\t{t.mt}\t{t.rrna}\n")


def write_share_rows(out, label, attr, tallies, categories, scope=None, available=True):
    """A count and a percentage row per category, for the *attr* ("rrna"/"mt") share.

    *scope* is appended to every row's qualifiers ("called cells"). Without
    *tallies*, or when not *available*, the counts are N/A rather than absent, so
    "not computed" is not read as 0.
    """
    for key, unit, pool, qualifier in categories:
        tally = None if tallies is None else tallies[key]
        total = NA if tally is None else tally.total
        part = getattr(tally, attr) if tally is not None and available else NA
        quals = [q for q in (qualifier, scope) if q]
        suffix = f" ({', '.join(quals)})" if quals else ""
        out.row(f"{label} {unit} in {pool}{suffix}", part)
        out.row(f"Percentage of {label} {unit} (of {', '.join([pool] + quals)})", frac(part, total))


def write_annotation_rows(out, args, ann, regions):
    out.text("GTF file", args.gtf)
    out.text("Biotype attribute used", ann["biotype_attr"] or "none")
    out.text("Feature type used", ann["feature_type"] or "none")
    out.text("rRNA biotypes detected", ann["rrna_biotypes"])
    out.text("rRNA reference contigs (added)", " ".join(regions.added_present))
    out.row("rRNA regions counted", regions.n_rrna_regions)
    out.text("MT Contig", " ".join(args.mt_contig))


def write_library_rows(out, c, regions):
    # Library-level rows carry no scope suffix: they describe every read in the BAM,
    # empty droplets, ambient RNA and barcodes outside the whitelist included
    library = c.library
    out.row("Mapped reads", library["mapped"].total)
    out.row("Unmapped reads", c.unmapped)
    out.row("Uniquely mapped reads", library["unique"].total)
    out.row("Multimapped reads (primary alignment)", library["mmpa"].total)
    out.row("Multimapped read alignments (all alignments)", library["mmaa"].total)

    # Mt_rRNA is an rRNA biotype, so reads on MT-RNR1/MT-RNR2 count here and as mtDNA
    # reads too. The "all alignments" rows show where multimapper alignments land, not
    # a fraction of reads.
    write_share_rows(out, "rRNA", "rrna", library, CATEGORIES, available=regions.n_rrna_regions > 0)
    # An mtDNA read is a read whose primary alignment is on a mitochondrial contig
    write_share_rows(out, "mtDNA", "mt", library, CATEGORIES)


def write_mt_origin_rows(out, c, regions):
    """mtDNA reads by origin.

    No read can be proven to come from DNA, and the mitochondrial genome is
    transcribed from both strands almost end to end, so these are indicators, not
    a contamination measurement: DNA-derived fragments land on both strands about
    equally and also cover untranscribed stretches, transcripts land sense to the
    annotated genes. A read where genes on opposite strands overlap (MT-ND5/MT-ND6)
    counts as both sense and antisense.
    """
    mt_total = mt_sense = mt_antisense = mt_outside = NA
    if regions.has_mt_genes:
        mt_total = c.library["mapped"].mt
        mt_outside = mt_total - c.mt_in_genes
        if regions.sense_on_read_strand is not None:
            mt_sense, mt_antisense = c.mt_sense, c.mt_antisense
    for origin, n in (("sense to", mt_sense), ("antisense to", mt_antisense), ("outside", mt_outside)):
        out.row(f"mtDNA reads {origin} mitochondrial genes (primary alignment)", n)
        out.row(f"Percentage of mtDNA reads {origin} mitochondrial genes (of mtDNA reads)",
                frac(n, mt_total))


def write_called_cell_rows(out, c, regions, cells):
    """The read-level metrics, restricted to the barcodes of the filtered matrix.

    Cell calling only: doublet removal and CellSweep come later. Without a list
    every row is N/A rather than absent, so "not computed" is not read as 0.
    """
    if cells is None:
        tallies, n_cells, n_found, n_mapped, n_unique = None, NA, NA, NA, NA
    else:
        tallies = c.in_cells
        n_cells, n_found = len(cells), len(c.cells_found)
        n_mapped, n_unique = tallies["mapped"].total, tallies["unique"].total
    out.row("Called cell barcodes", n_cells)
    out.row("Called cell barcodes found in BAM", n_found)
    out.row("Mapped reads (called cells)", n_mapped)
    out.row("Uniquely mapped reads (called cells)", n_unique)
    write_share_rows(out, "rRNA", "rrna", tallies, CELL_CATEGORIES, "called cells",
                     available=regions.n_rrna_regions > 0)
    write_share_rows(out, "mtDNA", "mt", tallies, CELL_CATEGORIES, "called cells")


def write_mt_rrna(args):
    """The mtDNA / rRNA metrics and the per-barcode reads; returns the called cells or None."""
    with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as bam:
        ann, regions = load_regions(args, bam)
        cells = load_cells(args.cell_barcodes)
        counts = scan_bam(bam, regions, cells)
    write_barcode_reads(args.out_barcode_reads, counts.per_barcode)

    if cells is not None and not counts.cells_found:
        warn(f"none of the {len(cells)} called-cell barcodes occurs as a CB tag in {args.bam};",
             "check that barcodes.tsv and the BAM use the same barcode format")

    with MetricsWriter(args.out_mt_rrna) as out:
        write_annotation_rows(out, args, ann, regions)
        write_library_rows(out, counts, regions)
        write_mt_origin_rows(out, counts, regions)
        write_called_cell_rows(out, counts, regions, cells)
    return cells


# ------------------------------------------------------------------
# Antisense, from STARsolo's CellReads.stats
# ------------------------------------------------------------------

ANTISENSE_COLUMNS = ("exonic", "intronic", "exonicAS", "intronicAS")


def sum_cellreads(path, cells):
    """``(library, in_cells, n_rows_cells)`` sums of ANTISENSE_COLUMNS in CellReads.stats.

    Every row, CBnotInPasslist included, adds to the library-level sums. None, with
    a warning, when the file lacks one of the columns.
    """
    library = Counter(dict.fromkeys(ANTISENSE_COLUMNS, 0))
    in_cells = Counter(dict.fromkeys(ANTISENSE_COLUMNS, 0))
    n_rows_cells = 0
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [col for col in ANTISENSE_COLUMNS if col not in (reader.fieldnames or [])]
        if missing:
            warn(f"{path} lacks column(s) {' '.join(missing)}; no antisense metrics")
            return None
        barcode_col = reader.fieldnames[0]
        for rec in reader:
            values = {col: int(rec[col]) for col in ANTISENSE_COLUMNS}
            library.update(values)
            if cells is not None and rec[barcode_col] in cells:
                n_rows_cells += 1
                in_cells.update(values)
    return library, in_cells, n_rows_cells


def write_antisense_rows(out, counts, suffix):
    """Sense and antisense reads and the antisense shares; *suffix* scopes the row names."""
    sense = counts["exonic"] + counts["intronic"]
    anti = counts["exonicAS"] + counts["intronicAS"]
    out.row(f"Reads sense to genes (exonic + intronic{suffix})", sense)
    out.row(f"Reads antisense to genes (exonicAS + intronicAS{suffix})", anti)
    out.row(f"Percentage of antisense reads (of reads assigned to genes{suffix})", frac(anti, sense + anti))
    for kind in ("exonic", "intronic"):
        out.row(f"Percentage of {kind} antisense reads (of {kind} reads{suffix})",
                frac(counts[kind + "AS"], counts[kind] + counts[kind + "AS"]))


def write_antisense(args, cells):
    """Antisense share of gene-overlapping reads, from CellReads.stats.

    STARsolo flags only reads with a unique genomic alignment, against the gene
    model of the feature the file belongs to (GeneFull_Ex50pAS), on the strand set
    by --soloStrand. The four flags are mutually exclusive. GeneFull_Ex50pAS leaves
    exonicAS reads out of the matrix; intronicAS reads may still be counted in it.
    """
    if args.strand not in ("Forward", "Reverse"):
        warn(f"strand '{args.strand or 'unset'}' is not Forward or Reverse; no antisense metrics")
        return
    if not args.cellreads or not os.path.isfile(args.cellreads):
        warn("no STARsolo CellReads.stats (soloCellReadStats off, or a bam_only run); no antisense metrics")
        return
    sums = sum_cellreads(args.cellreads, cells)
    if sums is None:
        return
    library, in_cells, n_rows_cells = sums

    with MetricsWriter(args.out_antisense) as out:
        out.text("Source", f"STARsolo {os.path.basename(args.cellreads)} (GeneFull_Ex50pAS), uniquely mapped reads")
        out.text("STARsolo strand", args.strand)
        write_antisense_rows(out, library, "")
        if cells is not None:
            out.row("Called cell barcodes found in CellReads.stats", n_rows_cells)
            write_antisense_rows(out, in_cells, ", called cells")


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True, help="STARsolo BAM (coordinate-sorted)")
    p.add_argument("--gtf", required=True, help="Annotation the reads were mapped with (merged with ref_gtf_addfeature)")
    p.add_argument("--mt-contig", nargs="*", default=[],
                   help="Mitochondrial contig names; also accepts one space-separated string")
    p.add_argument("--strand", default="", help="star_soloStrand: Forward, Reverse or Unstranded")
    p.add_argument("--rrna-gtf", default="", help="The added rRNA reference (ref_gtf_addfeature)")
    p.add_argument("--cell-barcodes", default="", help="barcodes.tsv of the filtered matrix; adds the called-cell rows")
    p.add_argument("--cellreads", default="", help="STARsolo GeneFull_Ex50pAS/CellReads.stats, for the antisense metrics")
    p.add_argument("--out-mt-rrna", required=True, help="Output: mtDNA / rRNA metrics")
    p.add_argument("--out-antisense", required=True, help="Output: antisense metrics (stranded runs only)")
    p.add_argument("--out-barcode-reads", required=True,
                   help="Output: mapped, mtDNA and rRNA reads per barcode (gzipped TSV)")
    p.add_argument("--threads", type=int, default=1, help="BGZF decompression threads")
    args = p.parse_args()

    # params.mt_contig may arrive as one quoted argument ("chrM M MT") or as words
    args.mt_contig = [c for arg in args.mt_contig for c in arg.split()]
    return args


def main():
    args = parse_args()
    cells = write_mt_rrna(args)
    write_antisense(args, cells)


if __name__ == "__main__":
    main()
