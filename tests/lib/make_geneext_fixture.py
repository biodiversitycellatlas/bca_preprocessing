#!/usr/bin/env python3
"""Write a synthetic GeneExt output directory for the dashboard tests.

GeneExt itself needs a genome, an annotation and an aligned BAM, none of which
exist on a developer machine, so the two files the dashboard reads back -- the
HTML report and the plain-text log -- are reconstructed here from the shapes
GeneExt emits.  The numbers are fixed and deliberately distinctive, so the check
can assert on them rather than on "something non-empty".

The report's ``ext_table`` is filled with one row per extended gene, which is
what makes GeneExt's own report large; the dashboard is expected to drop it.

Usage:
    make_geneext_fixture.py OUTDIR
"""

import json
import os
import random
import sys

# Numbers the check asserts on. Kept here so the two files agree with each
# other, except for the gene count, which GeneExt itself reports one lower in
# the log than in the report -- reproduced so the fallback path is told apart
# from the report path by its value alone.
N_GENES        = 31949
N_EXTENDED     = 13418
N_EXTENDED_LOG = 13417
MEDIAN_EXT     = 1383.0
MAX_EXT        = 2982.0
N_GENIC_PEAKS  = 89012
N_NOOV_PEAKS   = 167536
N_INITIAL_PEAKS = 244298
N_RETAINED_PEAKS = 128666
COV_PERCENTILE = 25
N_READS        = 396409600


def build_payload() -> dict:
    """The JSON object GeneExt embeds in its report as ``const D``."""
    rng = random.Random(0)
    return {
        "summary": {
            "n_genes": N_GENES, "n_extended": N_EXTENDED, "pct_extended": 42.0,
            "min_ext": 1.0, "median_ext": MEDIAN_EXT, "mean_ext": 1382.9,
            "max_ext": MAX_EXT,
            "n_genic_peaks": N_GENIC_PEAKS, "n_noov_peaks": N_NOOV_PEAKS,
            "n_orphan_peaks": 0, "orphan_warn_fraction": 0.1,
            "orphan_gene_fraction": 0.0, "orphan_warning": False,
            "cov_percentile": COV_PERCENTILE, "cov_threshold": "2.0",
            "n_reads": N_READS, "subsampled": False,
            "output_file": "geneext.gtf", "input_file": "genome.fixed.gtf",
            "genome_fixed": True, "log_genome_fix": True,
            "log_file": "geneext.gtf.GeneExt.log", "run_date": "2026-08-19 23:09",
            "run_args": "geneext.py -g ref.gtf -b merged_genome.bam -o geneext.gtf -j 4 -force",
        },
        "ext_hist": {
            "labels": [str(int(16 + i * 29.8)) for i in range(100)],
            "counts": [rng.randint(90, 500) for _ in range(100)],
        },
        "cov_hist": {
            "labels": [f"{-1.188 + i * 0.1145:.3f}" for i in range(50)],
            "counts_genic": [rng.randint(0, 8000) for _ in range(50)],
            "counts_noov": [rng.randint(0, 15000) for _ in range(50)],
            "log10_threshold": 0.301,
        },
        "mapping_stats": [],
        # The two blocks the dashboard is expected to leave behind
        "ext_table": [
            {"gene": f"Gene_{i:06d}", "peak": f"plus_peak_{i}", "ext": rng.randint(1, 2982)}
            for i in range(N_EXTENDED)
        ],
        "orphan_bed": "",
        "peak_flow": {
            "has_macs2_peaks": True,
            "initial_called": N_INITIAL_PEAKS,
            "passed_filtering": N_RETAINED_PEAKS,
            "assigned_to_genes": N_EXTENDED,
            "orphan_enabled": False, "orphan_count": 0,
            "filtered_file": "allpeaks_noov_fcov.bed",
        },
        "log_sections": ["Preflight checks", "Execution", "All done!"],
        "log_notes": [
            'Genome annotation warning: Could not find "gene" features in',
            "Retained 128666/167536 (76.8 %) intergenic peaks.",
        ],
        "fix_info": {
            "schema": "v1", "rerun_mode": False, "force_mode": True,
            "skipped_steps": [],
            "extension_param": {
                "name": "--maxdist", "mode": "auto", "user_value_bp": None,
                "effective_value_bp": 2982, "auto_quantile": 0.5,
            },
            "steps": {
                "mRNA_to_transcript": {"applied": False},
                "gene_features_added": {"applied": True, "n_genes_added": N_GENES},
                "clip_5prime": {"applied": False, "n_events": 0, "n_genes_clipped": 0},
            },
        },
    }


# The rich-rendered log GeneExt writes alongside the report. The box-drawing
# characters and the wrapped warning are kept: they are what a naive
# line-oriented parser trips over.
LOG_TEXT = f"""╭──────────────────╮
│ Preflight checks │
╰──────────────────╯
Genome annotation warning: Could not find "gene" features in
/path/to/reference/annotation.gtf! Trying to fix ...
╭───────────╮
│ Execution │
╰───────────╯
Running macs2 ... done
Filtering macs2 peaks ...
Retained {N_RETAINED_PEAKS}/{N_NOOV_PEAKS} (76.8 %) intergenic peaks.
done
Extending genes ... done
╭───────────╮
│ All done! │
╰───────────╯
Extended {N_EXTENDED_LOG}/{N_GENES} genes
Median extension length: {MEDIAN_EXT} bp
HTML report written to: geneext.gtf.Report.html
"""


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)

    report = (
        "<!DOCTYPE html>\n<html><head><title>GeneExt — geneext.gtf</title></head>\n"
        "<body>\n<main></main>\n<script>\n"
        "const D = " + json.dumps(build_payload()) + ";\n"
        "const S = D.summary;\n"
        "</script>\n</body></html>\n"
    )
    report_path = os.path.join(outdir, "geneext.gtf.Report.html")
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report)

    log_path = os.path.join(outdir, "geneext.gtf.GeneExt.log")
    with open(log_path, "w", encoding="utf-8") as fh:
        fh.write(LOG_TEXT)

    print(f"report: {report_path} ({os.path.getsize(report_path)} bytes)")
    print(f"log:    {log_path}")


if __name__ == "__main__":
    main()
