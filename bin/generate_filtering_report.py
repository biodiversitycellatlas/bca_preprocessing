#!/usr/bin/env python3
"""
Build ``filtering_report.html``: the interactive metacell filtering report.

Fills the template ``filtering_report.html`` with one JSON payload, assembled from

    *_mc_summary.json          run_metacells.py, one per sample
    *_mt_rrna_metrics.txt      CALC_READ_METRICS (BAM), STARsolo samples only
    *_antisense_metrics.txt    CALC_READ_METRICS (STARsolo CellReads.stats), STARsolo samples only
    CellReads.stats            STARsolo's per-barcode read counts, given as id=path pairs, for
                               the reads and sequencing saturation of every cell

the samplesheet's expected cells per sample (``--expected_cells``), shown beside the
called cells, and the default thresholds the report's sliders start at. The page does the
filtering itself and exports the selection ``apply_metacell_filter.py`` applies, so
nothing here decides what is kept.

The gene axis is the bulk of the payload, and samples counted against the same
annotation share it: every distinct gene axis becomes one "universe" (ids, names,
mito/rRNA flags, domains as indices into one run-wide dictionary), and each sample
carries only its own per-gene numbers aligned to it. Genes with no UMI in any sample
of a universe are left out of it; ``apply_metacell_filter.py`` excludes those by
rule rather than by list.

The header carries the BCA logo (``--logo``, assets/bca_logo.svg in the pipeline) in
the style of dashboard_report.html.

Standard library only, like generate_dashboard.py, so it runs without a container.
Every embedded block is escaped against ``</script>``, and the template is filled
in a single pass, so neither a gene name nor a domain can break out of its block or
be mistaken for a placeholder.
"""

import argparse
import base64
import csv
import datetime
import json
import logging
import os
import re
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("FilteringReport")

SUMMARY_SCHEMA = "bca_metacell_summary"
PLACEHOLDER_RE = re.compile(r"__([A-Z_]+)_PLACEHOLDER__")

# The suffixes apply_suffix in mapping_workflow.nf gives an analytical run's id
MAPPER_SUFFIXES = ("_geneext_subsampled_starsolo", "_geneext_alevinfry", "_geneext_starsolo",
                   "_subsampled_starsolo", "_starsolo", "_alevinfry")

# Rows of *_mt_rrna_metrics.txt the report shows: the called-cell value first, the
# whole-BAM value when no barcode list reached CALC_READ_METRICS. Both are percentages of
# mapped reads (primary alignment), so mtDNA and rRNA share one denominator. Rows of
# earlier pipeline versions are not read: their mtDNA value counted alignment
# records and their rRNA value had uniquely mapped reads as the denominator.
FC_KEYS = {
    "mt_reads_pct": (
        "Percentage of mtDNA reads (of mapped reads, primary alignment, called cells)",
        "Percentage of mtDNA reads (of mapped reads, primary alignment)",
    ),
    "rrna_reads_pct": (
        "Percentage of rRNA reads (of mapped reads, primary alignment, called cells)",
        "Percentage of rRNA reads (of mapped reads, primary alignment)",
    ),
}
# Antisense comes from STARsolo's CellReads.stats, summed over the called cells or, without
# a barcode list, over every barcode: the whole library rather than the BAM
ANTISENSE_KEYS = (
    "Percentage of antisense reads (of reads assigned to genes, called cells)",
    "Percentage of antisense reads (of reads assigned to genes)",
)

# Columns of CellReads.stats, as secondderiv_filter_matrices.py reads them: the reads counted
# into the matrix (unique-gene reads, in preference order) and their deduplicated UMIs
CELLREADS_READ_COLUMNS = ("countedU", "featureU")
CELLREADS_UMI_COLUMN = "nUMIunique"
# The share of a sample's cells that must be found in a CellReads.stats before its reads
# are used; below it the barcodes are not the same cells (another mapper's barcode form)
CELLREADS_MIN_MATCH = 0.9


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------

def parse_metrics(path):
    """A read metrics file (``key,value`` CSV rows) as a dict.

    Read as CSV: metric names carry commas ("of mapped reads, primary alignment")
    and are quoted. Files of earlier versions, unquoted, read the same way.
    """
    data = {}
    if not path or not os.path.exists(path):
        return data
    with open(path, newline="") as fh:
        for rec in csv.reader(fh):
            if len(rec) >= 2:
                data[rec[0].strip()] = ",".join(rec[1:]).strip()
    return data


def to_pct(value):
    """The scripts write fractions (0.0123); the report shows percentages."""
    try:
        return round(100.0 * float(value), 3)
    except (TypeError, ValueError):
        return None


def base_id(sample_id):
    for suffix in MAPPER_SUFFIXES:
        if sample_id.endswith(suffix):
            return sample_id[: -len(suffix)]
    return sample_id


def index_by_suffix(paths, suffix):
    """``{analytical_id: path}`` for files named ``<id><suffix>``, matched exactly."""
    out = {}
    for p in paths or []:
        name = os.path.basename(p)
        if name.endswith(suffix):
            out[name[: -len(suffix)]] = p
    return out


def pick_scoped(metrics, cells_key, bam_key, library_scope="whole BAM"):
    """``(percentage, scope)``: the called-cell value when there is one, else the
    library-level value (labelled *library_scope*), else ``(None, None)``."""
    for key, scope in ((cells_key, "called cells"), (bam_key, library_scope)):
        value = to_pct(metrics.get(key))
        if value is not None:
            return value, scope
    return None, None


def starsolo_candidates(sample_id):
    """
    The runs whose STARsolo outputs can stand in for this sample's: the sample itself,
    then the STARsolo run of the same sample. A GeneExt run takes the STARsolo run
    mapped against the same extended annotation first.
    """
    candidates = [sample_id]
    b = base_id(sample_id)
    if sample_id.endswith("_geneext_alevinfry"):
        candidates += [b + "_geneext_starsolo", b + "_geneext_subsampled_starsolo"]
    candidates += [b + "_starsolo", b + "_subsampled_starsolo"]
    return candidates


def featurecounts_for(sample_id, mt_files, antisense_files):
    """
    The sample's read metrics, or those of the STARsolo run of the same sample
    (they are read from STARsolo's BAM and CellReads.stats), labelled so.
    """
    for cid in starsolo_candidates(sample_id):
        if cid in mt_files or cid in antisense_files:
            mt = parse_metrics(mt_files.get(cid))
            anti = parse_metrics(antisense_files.get(cid))
            out = {
                "source_id": cid,
                "from_other_mapper": cid != sample_id,
            }
            out["antisense_pct"], out["antisense_scope"] = pick_scoped(anti, *ANTISENSE_KEYS, "whole library")
            for key, (cells_key, bam_key) in FC_KEYS.items():
                out[key], out[key.replace("_pct", "_scope")] = pick_scoped(mt, cells_key, bam_key)
            return out
    return None


def parse_pairs(pairs, option):
    """``{sample_id: value}`` from ``id=value`` pairs."""
    out = {}
    for pair in pairs or []:
        sample_id, sep, value = pair.rpartition("=")
        if not sep:
            raise ValueError(f"{option} takes id=value pairs, got {pair!r}")
        out[sample_id] = value
    return out


def parse_expected_cells(pairs):
    """``{sample_id: expected_cells}`` from ``id=n`` pairs: the samplesheet's numbers,
    for the samples that give one."""
    return {k: int(v) for k, v in parse_pairs(pairs, "--expected_cells").items()}


def parse_cellreads(path):
    """
    ``{barcode: (reads, umis)}`` from a STARsolo CellReads.stats: the reads counted into
    the matrix and their unique UMIs. Empty when the file lacks the columns.
    """
    out = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None) or []
        reads_col = next((header.index(c) for c in CELLREADS_READ_COLUMNS if c in header), None)
        if reads_col is None or CELLREADS_UMI_COLUMN not in header:
            logger.warning(f"{path} has none of {CELLREADS_READ_COLUMNS} or no {CELLREADS_UMI_COLUMN}; "
                           "no per-cell reads from it")
            return out
        umi_col = header.index(CELLREADS_UMI_COLUMN)
        for rec in reader:
            try:
                out[rec[0]] = (float(rec[reads_col]), float(rec[umi_col]))
            except (IndexError, ValueError):
                continue
    return out


def cell_reads_for(sample_id, barcodes, cellreads_files, cache):
    """
    ``(reads, saturation, source_id)`` per cell from the sample's CellReads.stats, or that
    of the STARsolo run of the same sample. A file is used only when it holds at least
    CELLREADS_MIN_MATCH of the sample's cells, as another mapper names barcodes its own way.
    Saturation is 1 - unique UMIs / counted reads; None where a cell has no counted read.
    *cache* keeps each parsed file for the other samples that borrow it.
    """
    for cid in starsolo_candidates(sample_id):
        if cid not in cellreads_files:
            continue
        if cid not in cache:
            cache[cid] = parse_cellreads(cellreads_files[cid])
        stats = cache[cid]
        found = [stats.get(b) for b in barcodes]
        n_found = sum(1 for f in found if f is not None)
        if not barcodes or n_found < CELLREADS_MIN_MATCH * len(barcodes):
            logger.info(f"{sample_id}: {cid}'s CellReads.stats holds {n_found} of {len(barcodes)} cells; not used")
            continue
        reads = [int(f[0]) if f else None for f in found]
        saturation = [round(1.0 - f[1] / f[0], 4) if f and f[0] > 0 else None for f in found]
        return reads, saturation, cid
    return None, None, None


def load_summaries(paths):
    summaries = []
    for p in paths or []:
        with open(p) as fh:
            s = json.load(fh)
        if s.get("schema") != SUMMARY_SCHEMA:
            logger.warning(f"Skipping {p}: not a metacell summary")
            continue
        summaries.append(s)
    summaries.sort(key=lambda s: s["sample"]["id"])
    return summaries


# --------------------------------------------------------------------------
# Payload
# --------------------------------------------------------------------------

def columnar(records, keys):
    return {k: [r.get(k) for r in records] for k in keys}


METACELL_KEYS = ["id", "pseudo", "n_cells", "umis", "median_umis", "mito_pct", "median_cell_mito_pct",
                 "rrna_pct", "n_doublets", "doublet_pct", "alpha_hat_mean", "x", "y", "markers"]
CELL_KEYS = ["barcodes", "group", "umis", "n_genes", "mito_umis", "rrna_umis", "doublet", "alpha"]


def build_universes(summaries):
    """
    One universe per distinct gene axis, each sample pointing at its own.

    Returns ``(universes, domains, per_sample)`` where ``per_sample[i]`` is
    ``(universe_index, kept_positions)`` for summaries[i], or None for a skipped one.
    """
    by_axis = {}
    groups = []
    for i, s in enumerate(summaries):
        genes = s.get("genes")
        if not genes:
            continue
        key = tuple(genes["ids"])
        if key not in by_axis:
            by_axis[key] = len(groups)
            groups.append([])
        groups[by_axis[key]].append(i)

    domain_index = {}
    domains = []
    universes = []
    per_sample = [None] * len(summaries)

    for members in groups:
        first = summaries[members[0]]["genes"]
        n = len(first["ids"])
        expressed = [False] * n
        for i in members:
            for j, v in enumerate(summaries[i]["genes"]["total_umi"]):
                if v:
                    expressed[j] = True
        keep = [j for j in range(n) if expressed[j]]

        pfam = []
        for j in keep:
            idx = []
            for d in first["pfam"][j]:
                if d not in domain_index:
                    domain_index[d] = len(domains)
                    domains.append(d)
                idx.append(domain_index[d])
            pfam.append(idx)

        universes.append({
            "ids": [first["ids"][j] for j in keep],
            "names": [first["names"][j] for j in keep],
            "mito": [first["mito"][j] for j in keep],
            "rrna": [first["rrna"][j] for j in keep],
            "pfam": pfam,
            "n_unexpressed": n - len(keep),
        })
        for i in members:
            per_sample[i] = (len(universes) - 1, keep)

    return universes, domains, per_sample


def cells_block(s, cellreads_files, cache):
    """
    The summary's per-cell arrays, with the reads and saturation of each cell from
    CellReads.stats where a matching one was given; None for a summary without them
    (written before per-cell arrays, or a skipped sample).
    """
    cells = s.get("cells")
    if not cells:
        return None
    out = {k: cells.get(k) for k in CELL_KEYS}
    out["reads"], out["saturation"], out["reads_source"] = cell_reads_for(
        s["sample"]["id"], cells["barcodes"], cellreads_files, cache)
    return out


def expression_block(expression, keep):
    """The summary's metacell profiles, their gene positions moved onto the universe's."""
    if not expression:
        return None
    position = {j: k for k, j in enumerate(keep)}
    rows = [(position[j], peak, q) for j, peak, q in zip(expression["genes"], expression["max"], expression["q"])
            if j in position]
    return {"mc": expression["mc"], "genes": [r[0] for r in rows],
            "max": [r[1] for r in rows], "q": [r[2] for r in rows]}


def build_payload(summaries, mt_files, antisense_files, defaults, run_info, expected_cells=None,
                  cellreads_files=None):
    universes, domains, per_sample = build_universes(summaries)
    expected_cells = expected_cells or {}
    cellreads_files = cellreads_files or {}
    cellreads_cache = {}

    samples = []
    for s, placement in zip(summaries, per_sample):
        sample = dict(s["sample"])
        sample["expected_cells"] = expected_cells.get(sample["id"])
        sample["featurecounts"] = featurecounts_for(sample["id"], mt_files, antisense_files)
        sample["metacells"] = columnar(s.get("metacells") or [], METACELL_KEYS)
        sample["cells"] = cells_block(s, cellreads_files, cellreads_cache)

        if placement is not None:
            u, keep = placement
            g = s["genes"]
            sample["universe"] = u
            sample["genes"] = {
                "total_umi": [g["total_umi"][j] for j in keep],
                "n_cells": [g["n_cells"][j] for j in keep],
                "ambient_hat": [g["ambient_hat"][j] for j in keep] if g.get("ambient_hat") else None,
                "tau": [g["tau"][j] for j in keep] if g.get("tau") else None,
            }
            sample["expression"] = expression_block(s.get("expression"), keep)
        else:
            sample["universe"] = None
            sample["genes"] = None
            sample["expression"] = None
        samples.append(sample)

    return {
        "run": run_info,
        "defaults": defaults,
        "domains": domains,
        "universes": universes,
        "samples": samples,
    }


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------

def safe_json(obj):
    """JSON that cannot close the <script> block it is embedded in."""
    return (json.dumps(obj, separators=(",", ":"), allow_nan=False)
            .replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
            .replace("\u2028", "\\u2028").replace("\u2029", "\\u2029"))


def logo_tag(path):
    """
    The BCA logo for the header, inlined so the report stays one portable file. It is
    the same image dashboard_report.html embeds, kept once in assets/bca_logo.svg.
    """
    if not path:
        return ""
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    return f'<img src="data:image/svg+xml;base64,{data}" class="logo" alt="Biodiversity Cell Atlas logo">'


def render(template, blocks):
    """
    Fill every ``__NAME_PLACEHOLDER__`` in one pass, so text substituted in is never
    scanned for placeholders again. An unknown placeholder is an error.
    """
    missing = set(PLACEHOLDER_RE.findall(template)) - set(blocks)
    if missing:
        raise KeyError(f"The template has placeholders with no value: {sorted(missing)}")
    return PLACEHOLDER_RE.sub(lambda m: blocks[m.group(1)], template)


def main():
    parser = argparse.ArgumentParser(description="Build the interactive metacell filtering report.")
    parser.add_argument("--template", required=True, help="filtering_report.html template")
    parser.add_argument("--summaries", nargs="*", default=[], help="*_mc_summary.json files")
    parser.add_argument("--mt_rrna_metrics", nargs="*", default=[], help="*_mt_rrna_metrics.txt files")
    parser.add_argument("--antisense_metrics", nargs="*", default=[], help="*_antisense_metrics.txt files")
    parser.add_argument("--expected_cells", nargs="*", default=[],
                        help="The samplesheet's expected cells as id=n pairs, for the samples that give one")
    parser.add_argument("--cellreads_stats", nargs="*", default=[],
                        help="STARsolo CellReads.stats files as id=path pairs, for per-cell reads and saturation")
    parser.add_argument("--max_mito_pct", type=float, default=None, help="Initial mito %% threshold")
    parser.add_argument("--max_doublet_pct", type=float, default=None, help="Initial doublet %% threshold")
    parser.add_argument("--max_alpha_hat", type=float, default=None, help="Initial ambient-fraction threshold")
    parser.add_argument("--min_gene_umis", type=float, default=None, help="Initial min total UMIs per gene")
    parser.add_argument("--pfam_presets", default="",
                        help="Domain regular expressions offered as presets, ';'-separated")
    parser.add_argument("--version", default="", help="Pipeline version")
    parser.add_argument("--logo", default=None,
                        help="Header logo (SVG), e.g. assets/bca_logo.svg; without it the header is text only")
    parser.add_argument("--output", default="filtering_report.html")
    args = parser.parse_args()

    summaries = load_summaries(args.summaries)
    if not summaries:
        logger.warning("No metacell summaries given; the report will say so")

    defaults = {
        "max_mito_pct": args.max_mito_pct,
        "max_doublet_pct": args.max_doublet_pct,
        "max_alpha_hat": args.max_alpha_hat,
        "min_gene_umis": args.min_gene_umis,
        "pfam_presets": [p for p in (s.strip() for s in args.pfam_presets.split(";")) if p],
    }
    run_info = {
        "pipeline_version": args.version,
        "generated": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    }

    payload = build_payload(
        summaries,
        index_by_suffix(args.mt_rrna_metrics, "_mt_rrna_metrics.txt"),
        index_by_suffix(args.antisense_metrics, "_antisense_metrics.txt"),
        defaults, run_info,
        parse_expected_cells(args.expected_cells),
        parse_pairs(args.cellreads_stats, "--cellreads_stats"),
    )

    with open(args.template, encoding="utf-8") as fh:
        template = fh.read()
    data = safe_json(payload)
    html = render(template, {"REPORT_DATA": data, "LOGO": logo_tag(args.logo)})

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(html)

    n_ok = sum(1 for s in payload["samples"] if s.get("status") == "ok")
    n_cells = sum(1 for s in payload["samples"] if s.get("cells"))
    n_reads = sum(1 for s in payload["samples"] if s.get("cells") and s["cells"]["reads"])
    logger.info(
        f"Wrote {args.output}: {len(payload['samples'])} samples ({n_ok} with metacells, {n_cells} with "
        f"per-cell data, {n_reads} with per-cell reads), {len(payload['universes'])} gene universe(s), "
        f"{len(payload['domains'])} domains; {len(data) / 1e6:.1f} MB of data"
    )


if __name__ == "__main__":
    main()
