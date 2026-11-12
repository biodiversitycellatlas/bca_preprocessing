#!/usr/bin/env python3
"""
One row per gene of the reference annotation, carrying what the metacell filtering
report filters on: whether the gene is mitochondrial, whether it is an rRNA, and the
protein domains a user-supplied annotation gives it, plus the name that annotation gives.

    gene_id  gene_name  chrom  biotype  is_mito  is_rrna  pfam  annot_name

``gene_name`` stays the GTF's, since the per-sample steps match name-keyed feature axes
against it and compare it with the matrix's own symbols; ``annot_name`` is the
annotation's (e.g. "CYP3A7" for a gene the GTF only calls "LOC116613690").

Built once per GTF and shared by every sample counted against it, so that the
per-sample steps match their feature axis against one table instead of each parsing
the GTF again.

Mitochondrial genes are the genes on ``params.mt_contig``, the same contigs
``CALC_READ_METRICS`` counts reads on; rRNA genes are those whose biotype matches
``params.grep_rrna``, by default the "rRNA" ``calculate_read_metrics.py`` looks for. Both are
decided per gene here, so the per-cell percentages downstream are UMI fractions of
annotated genes rather than the read fractions featureCounts reports.

The domain annotation is read in one of four formats, detected from its first lines
unless ``--annotation_format`` names one:

    pfamscan    raw pfam_scan.pl output: '#' comment block, whitespace-aligned
                columns, one row per hit; the key is <seq id>, the domain <hmm name>
    emapper     eggNOG-mapper's .emapper.annotations: '##' preamble and footer, a
                '#query' header; domains from PFAMs, the name from Preferred_name
                (unless --annotation_id_col / _pfam_col / _name_col name columns it has)
    table       a TSV whose header holds --annotation_id_col and --annotation_pfam_col
    headerless  a TSV without header: 'id <tab> domains', or 'id <tab> name <tab> domains'
                when any row has three columns (e.g. <sp>_long.pep.annotations.csv,
                <sp>_long.pep.pfamscan_archs.csv)

The key is whatever identifier the tool that produced the file used. eggNOG-mapper,
pfam_scan and InterProScan key on protein or transcript ids, so the key is resolved
back to a gene through the GTF, in order:

    gene_id -> transcript_id -> protein_id -> gene_name
    -> the same again with a trailing version (".1") or TransDecoder (".p1") suffix removed

Domains are treated as opaque tokens -- Pfam accessions (``PF00125``, version
stripped) or names (``Histone``, as eggNOG-mapper's ``PFAMs`` column has them),
separated by ``,`` ``;`` ``|`` ``/`` or whitespace (a domain architecture such as
``SAM_2/EGF/`` included) -- and a gene keeps the union over every row that resolved to
it. Its name is the first non-empty one among those rows.
"""

import argparse
import json
import logging
import re
import sys
from collections import OrderedDict, defaultdict

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("BuildGeneTable")

GENE_TABLE_COLUMNS = ["gene_id", "gene_name", "chrom", "biotype", "is_mito", "is_rrna", "pfam", "annot_name"]
ANNOTATION_FORMATS = ("auto", "table", "emapper", "pfamscan", "headerless")

_ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')
_PFAM_SPLIT_RE = re.compile(r"[,;|/\s]+")
_PFAM_VERSION_RE = re.compile(r"^(PF\d{5})\.\d+$")
_EMPTY_TOKENS = {"", "-", "na", "nan", "none", "null"}

# Biotype attributes, in the order GENCODE, Ensembl and RefSeq-style GTFs spell them
_GENE_BIOTYPE_KEYS = ("gene_biotype", "gene_type")
_TX_BIOTYPE_KEYS = ("transcript_biotype", "transcript_type")


def parse_attributes(field):
    """A GTF attribute column as a dict, first occurrence of each key winning."""
    attrs = {}
    for key, value in _ATTR_RE.findall(field):
        attrs.setdefault(key, value)
    return attrs


def parse_gtf(gtf_path):
    """
    Genes and the identifiers that point back to them.

    Returns ``(genes, aliases)``: ``genes`` maps gene_id to name, chrom and biotype,
    in file order; ``aliases`` maps transcript and protein ids to their gene_id.
    Works on GTFs without ``gene`` rows, which some non-model annotations omit.
    """
    genes = OrderedDict()
    aliases = {}

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            chrom, ftype, attrs = cols[0], cols[2], parse_attributes(cols[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue

            entry = genes.get(gene_id)
            if entry is None:
                entry = genes[gene_id] = {"gene_name": "", "chrom": chrom, "biotype": "", "ftypes": set()}
            entry["ftypes"].add(ftype)

            if not entry["gene_name"] and attrs.get("gene_name"):
                entry["gene_name"] = attrs["gene_name"]

            biotype = next((attrs[k] for k in _GENE_BIOTYPE_KEYS if k in attrs), "")
            if not biotype:
                biotype = next((attrs[k] for k in _TX_BIOTYPE_KEYS if k in attrs), "")
            if biotype and (ftype == "gene" or not entry["biotype"]):
                entry["biotype"] = biotype

            for key in ("transcript_id", "protein_id"):
                if attrs.get(key):
                    aliases.setdefault(attrs[key], gene_id)

    return genes, aliases


def flag_mito(genes, mt_contigs):
    contigs = set(mt_contigs)
    return {gid: g["chrom"] in contigs for gid, g in genes.items()}


def flag_rrna(genes, rrna_pattern):
    """
    rRNA by biotype, matched the way ``calculate_read_metrics.py`` looks for it, or by
    an explicit ``rRNA`` feature type.
    """
    pattern = re.compile(rrna_pattern) if rrna_pattern else None
    return {
        gid: bool((pattern and pattern.search(g["biotype"])) or "rRNA" in g["ftypes"])
        for gid, g in genes.items()
    }


def split_domains(cell):
    """The domain tokens of one annotation cell, versions stripped, empties dropped."""
    tokens = []
    for token in _PFAM_SPLIT_RE.split(str(cell)):
        token = token.strip()
        if token.lower() in _EMPTY_TOKENS:
            continue
        tokens.append(_PFAM_VERSION_RE.sub(r"\1", token))
    return tokens


def _key_variants(key):
    """The key as given, then with a version or TransDecoder suffix removed."""
    yield key
    stripped = re.sub(r"\.p\d+$", "", key)
    if stripped != key:
        yield stripped
    unversioned = re.sub(r"\.\d+$", "", stripped)
    if unversioned != stripped:
        yield unversioned


def resolver(genes, aliases):
    """
    A function mapping an annotation key to a gene_id, or None.

    Gene names are only used when they are unambiguous: a name shared by two genes
    would otherwise hand one of them the other's domains.
    """
    name_counts = defaultdict(int)
    for g in genes.values():
        if g["gene_name"]:
            name_counts[g["gene_name"]] += 1
    by_name = {g["gene_name"]: gid for gid, g in genes.items()
               if g["gene_name"] and name_counts[g["gene_name"]] == 1}

    def resolve(key):
        for variant in _key_variants(key):
            if variant in genes:
                return variant
            if variant in aliases:
                return aliases[variant]
            if variant in by_name:
                return by_name[variant]
        return None

    return resolve


def _header_columns(line):
    """A header line's column names; eggNOG-mapper writes its first one as '#query'."""
    return [c.strip().lstrip("#") for c in line.rstrip("\r\n").split("\t")]


def _is_pfamscan_comment(line):
    return line.startswith("#") and ("pfam_scan" in line or "<seq id>" in line)


def detect_format(annotation_path, id_col):
    """
    The annotation's format, from its first line that is neither blank nor a '##'
    comment. The launch check in utils_nfcore_bca_pipeline/main.nf applies the same rules.
    """
    with open(annotation_path, errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            if _is_pfamscan_comment(line):
                return "pfamscan"
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                return "emapper"
            return "table" if id_col in _header_columns(line) else "headerless"
    return "headerless"


def _column(header, path, wanted, fallback=None, required=True):
    """Index of the first of ``wanted``, ``fallback`` the header has."""
    for col in (wanted, fallback):
        if col and col in header:
            return header.index(col)
    if not required:
        return None
    raise ValueError(f"Column '{wanted}' not found in {path}. Available columns: {header}")


def _read_header_table(path, fmt, id_col, pfam_col, name_col):
    """
    Rows of a table with a header line; '##' lines (eggNOG-mapper's preamble and
    footer) are skipped wherever they are.

    Raises when a configured column is missing: a misnamed column would otherwise
    read as "no gene has a domain" and the report would just hide its PFAM panel.
    """
    emapper = fmt == "emapper"
    with open(path, errors="replace") as fh:
        idx = None
        for line in fh:
            if not line.strip() or line.startswith("##"):
                continue
            if idx is None:
                header = _header_columns(line)
                idx = (
                    _column(header, path, id_col, "query" if emapper else None),
                    _column(header, path, pfam_col, "PFAMs" if emapper else None),
                    _column(header, path, name_col) if name_col
                    else _column(header, path, "Preferred_name", required=False) if emapper else None,
                )
                continue
            fields = line.rstrip("\r\n").split("\t")
            # (key, name, domains), short rows padded
            yield tuple(fields[i] if i is not None and i < len(fields) else "" for i in (idx[0], idx[2], idx[1]))


def _read_pfamscan(path):
    """One row per pfam_scan.pl hit: <seq id> and <hmm name> (columns 1 and 7)."""
    with open(path, errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 7:
                raise ValueError(f"{path}, line {n}: expected pfam_scan.pl's 15 columns, found {len(fields)}")
            yield fields[0], "", fields[6]


def _read_headerless(path, id_col):
    """
    Rows of a header-less TSV: 'id <tab> domains', or 'id <tab> name <tab> domains'
    when any row has three columns -- decided for the whole file, so that a row whose
    trailing domain column is empty does not shift its name into the domain column.
    """
    rows = []
    with open(path, errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if any(f.strip() for f in fields[3:]):
                raise ValueError(
                    f"{path}, line {n}: {len(fields)} columns and no header with column '{id_col}'. "
                    "A file without header has 'id, domains' or 'id, name, domains'; otherwise set "
                    "annotation_id_col and annotation_pfam_col to columns of its header."
                )
            rows.append(fields[:3])
    width = max((len(r) for r in rows), default=1)
    for r in rows:
        r = r + [""] * (3 - len(r))
        yield (r[0], r[1], r[2]) if width == 3 else (r[0], "", r[1])


def iter_annotation_rows(path, fmt, id_col, pfam_col, name_col=None):
    """``(key, name, domain cell)`` per row of the annotation, in its ``fmt``."""
    if fmt == "pfamscan":
        return _read_pfamscan(path)
    if fmt == "headerless":
        return _read_headerless(path, id_col)
    return _read_header_table(path, fmt, id_col, pfam_col, name_col)


def load_domains(annotation_path, fmt, id_col, pfam_col, name_col, resolve):
    """Domains and name per gene_id, plus how the annotation rows resolved."""
    domains = defaultdict(set)
    names = {}
    n_rows = n_resolved = n_name_conflicts = 0
    unresolved_examples = []

    for key, name, cell in iter_annotation_rows(annotation_path, fmt, id_col, pfam_col, name_col):
        key = key.strip()
        if not key:
            continue
        n_rows += 1
        gene_id = resolve(key)
        if gene_id is None:
            if len(unresolved_examples) < 5:
                unresolved_examples.append(key)
            continue
        n_resolved += 1
        domains[gene_id].update(split_domains(cell))

        name = name.strip()
        if name.lower() not in _EMPTY_TOKENS:
            if gene_id not in names:
                names[gene_id] = name
            elif names[gene_id] != name:
                n_name_conflicts += 1

    stats = {
        "annotation_format": fmt,
        "annotation_rows": n_rows,
        "annotation_rows_resolved": n_resolved,
        "unresolved_examples": unresolved_examples,
        "annotation_name_conflicts": n_name_conflicts,
    }
    return domains, names, stats


def build_table(genes, is_mito, is_rrna, domains, names=None):
    names = names or {}
    rows = []
    for gid, g in genes.items():
        rows.append({
            "gene_id": gid,
            "gene_name": g["gene_name"],
            "chrom": g["chrom"],
            "biotype": g["biotype"],
            "is_mito": int(is_mito[gid]),
            "is_rrna": int(is_rrna[gid]),
            "pfam": ",".join(sorted(domains.get(gid, ()))),
            "annot_name": names.get(gid, ""),
        })
    return pd.DataFrame(rows, columns=GENE_TABLE_COLUMNS)


def main():
    parser = argparse.ArgumentParser(
        description="Build the per-gene table (mito / rRNA flags, protein domains) the metacell filtering uses."
    )
    parser.add_argument("--gtf", required=True, help="Reference GTF the matrices were counted against")
    parser.add_argument("--mt_contig", nargs="*", default=[],
                        help="Mitochondrial contig names (params.mt_contig, space-separated)")
    parser.add_argument("--rrna_pattern", default="rRNA",
                        help="Regular expression matched against the gene biotype (params.grep_rrna)")
    parser.add_argument("--annotation", default=None,
                        help="Optional protein domain annotation per gene/transcript/protein id "
                             "(pfam_scan.pl output, .emapper.annotations, a TSV with or without header)")
    parser.add_argument("--annotation_format", default="auto", choices=ANNOTATION_FORMATS,
                        help="Format of --annotation; 'auto' detects it from the first lines")
    parser.add_argument("--annotation_id_col", default="gene_id",
                        help="Identifier column of a --annotation with header ('query' is the emapper fallback)")
    parser.add_argument("--annotation_pfam_col", default="pfam",
                        help="Domain column of a --annotation with header ('PFAMs' is the emapper fallback)")
    parser.add_argument("--annotation_name_col", default=None,
                        help="Gene name column of a --annotation with header (emapper: Preferred_name)")
    parser.add_argument("--output", default="gene_table.tsv", help="Output TSV")
    parser.add_argument("--stats", default="gene_table_stats.json", help="Output JSON with match statistics")
    args = parser.parse_args()

    genes, aliases = parse_gtf(args.gtf)
    if not genes:
        raise SystemExit(f"No gene_id attributes found in {args.gtf}")
    logger.info(f"Read {len(genes)} genes and {len(aliases)} transcript/protein ids from {args.gtf}")

    is_mito = flag_mito(genes, args.mt_contig)
    is_rrna = flag_rrna(genes, args.rrna_pattern)

    warnings = []
    n_mito = sum(is_mito.values())
    if not args.mt_contig:
        warnings.append("No mitochondrial contig given (mt_contig is empty); mito percentages will be 0.")
    elif n_mito == 0:
        warnings.append(
            f"No genes found on the mitochondrial contig(s) {args.mt_contig}; mito percentages will be 0. "
            "Check mt_contig against the GTF's sequence names."
        )

    domains, names, ann_stats = {}, {}, {"annotation_format": None}
    if args.annotation:
        fmt = args.annotation_format
        if fmt == "auto":
            fmt = detect_format(args.annotation, args.annotation_id_col)
            logger.info(f"Annotation format detected as '{fmt}'")
        domains, names, ann_stats = load_domains(
            args.annotation, fmt, args.annotation_id_col, args.annotation_pfam_col,
            args.annotation_name_col, resolver(genes, aliases)
        )
        rate = ann_stats["annotation_rows_resolved"] / max(ann_stats["annotation_rows"], 1)
        logger.info(
            f"Resolved {ann_stats['annotation_rows_resolved']} / {ann_stats['annotation_rows']} annotation rows "
            f"({rate:.1%}) to {len(domains)} genes, {len(names)} of them named"
        )
        if ann_stats["annotation_name_conflicts"]:
            logger.info(
                f"{ann_stats['annotation_name_conflicts']} rows named their gene differently from an earlier "
                "row; the first name was kept"
            )
        if rate < 0.5:
            warnings.append(
                f"Only {rate:.1%} of the domain annotation's rows matched a gene, transcript or protein id "
                f"of the GTF (e.g. {ann_stats['unresolved_examples']}). Check annotation_id_col."
            )

    for w in warnings:
        logger.warning(w)

    table = build_table(genes, is_mito, is_rrna, domains, names)
    table.to_csv(args.output, sep="\t", index=False)

    stats = {
        "gtf": args.gtf,
        "n_genes": len(table),
        "n_mito": int(table["is_mito"].sum()),
        "n_rrna": int(table["is_rrna"].sum()),
        "mt_contig": args.mt_contig,
        "rrna_pattern": args.rrna_pattern,
        "annotation": args.annotation,
        "n_genes_with_domains": int((table["pfam"] != "").sum()),
        "n_genes_with_names": int((table["annot_name"] != "").sum()),
        "n_domains": len({d for ds in domains.values() for d in ds}),
        **ann_stats,
        "warnings": warnings,
    }
    with open(args.stats, "w") as fh:
        json.dump(stats, fh, indent=2)

    logger.info(
        f"Wrote {args.output}: {stats['n_genes']} genes, {stats['n_mito']} mitochondrial, "
        f"{stats['n_rrna']} rRNA, {stats['n_genes_with_domains']} with domains, "
        f"{stats['n_genes_with_names']} with an annotation name"
    )


if __name__ == "__main__":
    main()
