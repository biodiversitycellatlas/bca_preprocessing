#!/usr/bin/env python3
"""
One row per gene of the reference annotation, carrying what the metacell filtering
report filters on: whether the gene is mitochondrial, whether it is an rRNA, and the
protein domains a user-supplied annotation gives it.

    gene_id  gene_name  chrom  biotype  is_mito  is_rrna  pfam

Built once per GTF and shared by every sample counted against it, so that the
per-sample steps match their feature axis against one table instead of each parsing
the GTF again.

Mitochondrial genes are the genes on ``params.mt_contig``, the same contigs
``CALC_READ_METRICS`` counts reads on; rRNA genes are those whose biotype matches
``params.grep_rrna``, by default the "rRNA" ``calculate_read_metrics.py`` looks for. Both are
decided per gene here, so the per-cell percentages downstream are UMI fractions of
annotated genes rather than the read fractions featureCounts reports.

The domain annotation is a TSV with a header, keyed by whatever identifier the tool
that produced it used. eggNOG-mapper and InterProScan key on protein or transcript
ids, so the key is resolved back to a gene through the GTF, in order:

    gene_id -> transcript_id -> protein_id -> gene_name
    -> the same again with a trailing version (".1") or TransDecoder (".p1") suffix removed

Domains are treated as opaque tokens -- Pfam accessions (``PF00125``, version
stripped) or names (``Histone``, as eggNOG-mapper's ``PFAMs`` column has them) -- and
a gene keeps the union over every row that resolved to it.
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

GENE_TABLE_COLUMNS = ["gene_id", "gene_name", "chrom", "biotype", "is_mito", "is_rrna", "pfam"]

_ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')
_PFAM_SPLIT_RE = re.compile(r"[,;|\s]+")
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


def load_domains(annotation_path, id_col, pfam_col, resolve):
    """
    Domains per gene_id, plus how the annotation rows resolved.

    Raises when a configured column is missing: a misnamed column would otherwise
    read as "no gene has a domain" and the report would just hide its PFAM panel.
    """
    table = pd.read_csv(annotation_path, sep="\t", dtype=str, comment=None).fillna("")
    # eggNOG-mapper writes its header as '#query'
    table.columns = [c.lstrip("#") for c in table.columns]

    for col in (id_col, pfam_col):
        if col not in table.columns:
            raise ValueError(
                f"Column '{col}' not found in {annotation_path}. Available columns: {list(table.columns)}"
            )

    domains = defaultdict(set)
    n_rows = n_resolved = 0
    unresolved_examples = []

    for key, cell in zip(table[id_col], table[pfam_col]):
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

    stats = {
        "annotation_rows": n_rows,
        "annotation_rows_resolved": n_resolved,
        "unresolved_examples": unresolved_examples,
    }
    return domains, stats


def build_table(genes, is_mito, is_rrna, domains):
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
                        help="Optional TSV with a header giving protein domains per gene/transcript/protein id")
    parser.add_argument("--annotation_id_col", default="gene_id", help="Identifier column of --annotation")
    parser.add_argument("--annotation_pfam_col", default="pfam", help="Domain column of --annotation")
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

    domains, ann_stats = {}, {}
    if args.annotation:
        domains, ann_stats = load_domains(
            args.annotation, args.annotation_id_col, args.annotation_pfam_col, resolver(genes, aliases)
        )
        rate = ann_stats["annotation_rows_resolved"] / max(ann_stats["annotation_rows"], 1)
        logger.info(
            f"Resolved {ann_stats['annotation_rows_resolved']} / {ann_stats['annotation_rows']} annotation rows "
            f"({rate:.1%}) to {len(domains)} genes"
        )
        if rate < 0.5:
            warnings.append(
                f"Only {rate:.1%} of the domain annotation's rows matched a gene, transcript or protein id "
                f"of the GTF (e.g. {ann_stats['unresolved_examples']}). Check annotation_id_col."
            )

    for w in warnings:
        logger.warning(w)

    table = build_table(genes, is_mito, is_rrna, domains)
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
        "n_domains": len({d for ds in domains.values() for d in ds}),
        **ann_stats,
        "warnings": warnings,
    }
    with open(args.stats, "w") as fh:
        json.dump(stats, fh, indent=2)

    logger.info(
        f"Wrote {args.output}: {stats['n_genes']} genes, {stats['n_mito']} mitochondrial, "
        f"{stats['n_rrna']} rRNA, {stats['n_genes_with_domains']} with domains"
    )


if __name__ == "__main__":
    main()
