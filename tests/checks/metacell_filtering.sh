#!/usr/bin/env bash
# description: Check the metacell QC, the filtering report and the selection it exports, end to end
#
# The metacell filtering is a loop through a browser: run_metacells.py summarises each
# metacell, filtering_report.html turns the user's thresholds and lists into a
# selection, and apply_metacell_filter.py writes the matrices that selection asks for.
# Every link fails silently. A mito percentage computed after the mitochondrial genes
# were excluded reads 0 everywhere; an annotation lifted by position lands on the
# wrong cells; a selection applied to a regrouped sample keeps a plausible set of the
# wrong cells. None of it raises.
#
# tests/lib/make_metacell_fixture.py therefore plants the assignment (so Metacell2
# itself is not needed) and every number the report shows, and the cases assert on
# those values. Where a JavaScript runtime is available, the report's own script is
# run by tests/lib/report_harness.js to produce the selection, so the apply step is
# checked against what the page really exports.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

GENE_TABLE_TOOL="$PROJECT_ROOT/bin/build_gene_table.py"
SUMMARY_TOOL="$PROJECT_ROOT/bin/metacell_qc_summary.py"
RUN_TOOL="$PROJECT_ROOT/bin/run_metacells.py"
REPORT_TOOL="$PROJECT_ROOT/bin/generate_filtering_report.py"
APPLY_TOOL="$PROJECT_ROOT/bin/apply_metacell_filter.py"
TEMPLATE="$PROJECT_ROOT/bin/filtering_report.html"
FIXTURE_GEN="$TESTS_DIR/lib/make_metacell_fixture.py"
HARNESS="$TESTS_DIR/lib/report_harness.js"

PYTHON="${BCA_PYTHON:-}"
NODE="${BCA_NODE:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh metacell_filtering [-- OPTIONS]
       tests/checks/metacell_filtering.sh [OPTIONS]

Validate the metacell summaries, the filtering report and the apply step.

Options:
  --python PATH   Python interpreter to use (default: \$BCA_PYTHON, else
                  python3, else python). Needs numpy, pandas, scipy, anndata.
  --node PATH     JavaScript runtime for the report's script (default: \$BCA_NODE,
                  else node). Without one the report_js cases are skipped and the
                  selection is built in Python instead.
  --keep          Keep the generated fixtures and outputs for inspection.
  -h, --help      Show this message.

Cases:
  help                  the tools start and print usage
  gene_table            mito by contig, rRNA by biotype, domains via transcript and
                        protein ids with versions stripped
  match_genes           gene id, gene name and alevin-fry's cut id all match
  summary_qc            pooled mito %, doublet % and the lifted alpha_hat per metacell
  summary_modes         doublet % is null, not 0, when calls are absent or removed
  cells_h5ad            the per-cell annotations and CellSweep's layer land by barcode
  fingerprint           the fingerprint is stable and stored with the cells
  report_escape         every JSON block parses; no gene name can close the script
  report_payload        one shared gene universe; featureCounts as percentages,
                        called cells preferred, borrowed from STARsolo for alevin-fry
  report                featureCounts falls back to whole BAM, ignores legacy rows,
                        and pairs GeneExt runs
  report_js             the page's own filtering yields the expected selection and
                        survives a download -> load round trip
  apply_keep            exactly the kept metacells' cells, counts unchanged
  apply_genes           listed and unexpressed genes dropped, or flagged
  apply_guards          a stale fingerprint, unknown metacell or absent sample fails
  mc2_smoke             Metacell2 groups the fixture, deterministically (needs metacells)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --python)  PYTHON="${2:?--python needs a value}"; shift 2 ;;
        --node)    NODE="${2:?--node needs a value}"; shift 2 ;;
        --keep)    KEEP=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

# --------------------------------------------------------------------------
# Pre-flight
# --------------------------------------------------------------------------

log_header "metacell_filtering"

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"
    elif have_cmd python; then PYTHON="python"
    fi
fi

if [[ -z "$PYTHON" ]] || ! "$PYTHON" -c 'import sys; sys.exit(0)' >/dev/null 2>&1; then
    record SKIP "python" "no usable python interpreter found"
    finish_check
    exit 0
fi

if ! "$PYTHON" -c 'import numpy, pandas, scipy, anndata, h5py' >/dev/null 2>&1; then
    record SKIP "python" "numpy, pandas, scipy, anndata and h5py are needed"
    finish_check
    exit 0
fi

if [[ -z "$NODE" ]] && have_cmd node; then NODE="node"; fi
HAVE_NODE=0
[[ -n "$NODE" ]] && "$NODE" -e 'process.exit(0)' >/dev/null 2>&1 && HAVE_NODE=1

WORKDIR="$BCA_TEST_LOGDIR/metacell_filtering"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# assert NAME PYTHON_EXPR DETAIL [PATH...]
#
# As in annotated_h5ad.sh: PYTHON_EXPR yields (ok, message), with the fixture module
# and bin/ importable, and paths passed as A[0], A[1], ... so a Windows interpreter
# under Git Bash receives them converted.
assert() {
    local name="$1" expr="$2" detail="${3:-}"; shift 3 || shift $#
    local out
    if out="$("$PYTHON" -c "
import sys, json
LIB, BIN, W, A = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
sys.path.insert(0, LIB)
sys.path.insert(0, BIN)
import numpy as np, pandas as pd
import metacell_qc_summary as qcs
import make_metacell_fixture as F

def load(path):
    with open(path, encoding='utf-8') as fh:
        return json.load(fh)

def read_h5ad(path):
    import anndata as ad
    return ad.read_h5ad(path)

def dense(m):
    return np.asarray(m.todense()) if hasattr(m, 'todense') else np.asarray(m)

def mc(summary, group):
    return next(r for r in summary['metacells'] if r['id'] == group)

def close(a, b, tol=1e-3):
    return a is not None and b is not None and abs(a - b) <= tol

ok, msg = ($expr)
print(('PASS' if ok else 'FAIL') + '\t' + str(msg))
" "$TESTS_DIR/lib" "$PROJECT_ROOT/bin" "$WORKDIR" "$@" 2>&1)"; then
        local status="${out%%$'\t'*}" msg="${out#*$'\t'}"
        if [[ "$status" == "PASS" ]]; then
            record PASS "$name" "${detail:-$msg}"
        else
            record FAIL "$name" "$msg"
        fi
    else
        record FAIL "$name" "assertion crashed: $(printf '%s' "$out" | tail -3 | tr '\n' ' ')"
    fi
}

# step NAME LOG COMMAND... -- run a tool; a failure ends the check, since every later
# case reads what it wrote.
step() {
    local name="$1" logname="$2"; shift 2
    if run_logged "$BCA_TEST_LOGDIR/metacell_filtering_${logname}.log" "$@"; then
        return 0
    fi
    record FAIL "$name" "see metacell_filtering_${logname}.log"
    tail_log "$BCA_TEST_LOGDIR/metacell_filtering_${logname}.log"
    finish_check
    exit 1
}

# --------------------------------------------------------------------------
# Case: help
# --------------------------------------------------------------------------

HELP_OK=1
for tool in "$GENE_TABLE_TOOL" "$SUMMARY_TOOL" "$RUN_TOOL" "$REPORT_TOOL" "$APPLY_TOOL"; do
    name="$(basename "$tool")"
    if ! run_logged "$BCA_TEST_LOGDIR/metacell_filtering_help_$(slugify "$name").log" "$PYTHON" "$tool" --help; then
        record FAIL "help.$name" "--help failed"
        HELP_OK=0
    fi
done
if [[ "$HELP_OK" -eq 1 ]]; then
    record PASS "help" "5 tools print usage"
else
    finish_check
    exit 1
fi

# --------------------------------------------------------------------------
# Fixture and the gene table
# --------------------------------------------------------------------------

FIX="$WORKDIR/fixture"
step "fixture" fixture "$PYTHON" "$FIXTURE_GEN" "$FIX"
record PASS "fixture" "planted assignment, raw object, GTF and domain annotation written"

step "gene_table.run" gene_table "$PYTHON" "$GENE_TABLE_TOOL" \
    --gtf "$FIX/genes.gtf" --mt_contig chrM MT --rrna_pattern rRNA \
    --annotation "$FIX/pfam.tsv" --annotation_id_col query --annotation_pfam_col PFAMs \
    --output "$WORKDIR/gene_table.tsv" --stats "$WORKDIR/gene_table_stats.json"

assert "gene_table.flags" \
    "((lambda t: (sorted(t.gene_id[t.is_mito]) == sorted(F.MITO_GENES) and list(t.gene_id[t.is_rrna]) == F.RRNA_GENES
                  and list(t.gene_id) == F.GENE_IDS,
                  f'mito {list(t.gene_id[t.is_mito])}, rRNA {list(t.gene_id[t.is_rrna])}'))(qcs.load_gene_table(A[0])))" \
    "chrM genes are mitochondrial, the rRNA biotype is rRNA, GTF order kept" \
    "$WORKDIR/gene_table.tsv"

assert "gene_table.domains" \
    "((lambda t: ({g: d for g, d in zip(t.gene_id, t.pfam) if d} == F.EXPECTED_DOMAINS,
                  f\"domains {({g: d for g, d in zip(t.gene_id, t.pfam) if d})}\"))(qcs.load_gene_table(A[0])))" \
    "a versioned transcript key and a protein key resolve; PF00125.27 -> PF00125" \
    "$WORKDIR/gene_table.tsv"

assert "gene_table.stats" \
    "((lambda s: (s['annotation_rows'] == 4 and s['annotation_rows_resolved'] == 3
                  and s['unresolved_examples'] == ['UNKNOWN_KEY'] and s['n_mito'] == 2,
                  f\"{s['annotation_rows_resolved']} / {s['annotation_rows']} rows resolved\"))(load(A[0])))" \
    "the unresolvable key is counted and named" \
    "$WORKDIR/gene_table_stats.json"

# alevin-fry's ids are the GTF id cut at its first '-'; a name-only axis matches by name
assert "match_genes" \
    "((lambda r: (list(r[0]) == [0, 1, 4, -1] and r[1]['gene_id'] == 1 and r[1]['gene_name'] == 1
                  and r[1]['cut_id'] == 1 and r[1]['unmatched'] == 1,
                  f'rows {list(r[0])}, {r[1]}'))(
        qcs.match_genes(['GENE001', 'Gapdh', 'LOC5', 'NOT_A_GENE'], qcs.load_gene_table(A[0]))))" \
    "gene id, then gene name, then the cut id" \
    "$WORKDIR/gene_table.tsv"

# --------------------------------------------------------------------------
# run_metacells.py on the planted assignment
# --------------------------------------------------------------------------

OUT="$WORKDIR/mc"
mkdir -p "$OUT"
for sid in S1_starsolo S1_alevinfry; do
    step "run_metacells.$sid" "run_$sid" "$PYTHON" "$RUN_TOOL" \
        --cells_h5ad "$FIX/cells.h5ad" --ambient_h5ad "$FIX/raw.h5ad" \
        --gene_table "$WORKDIR/gene_table.tsv" --sample_id "$sid" --mapping_method fixture \
        --prefix "$OUT/$sid" --metacell_obs_col mc_truth --compression none
done
step "run_metacells.removed" run_removed "$PYTHON" "$RUN_TOOL" \
    --cells_h5ad "$FIX/cells.h5ad" --gene_table "$WORKDIR/gene_table.tsv" --sample_id S1_removed \
    --prefix "$OUT/S1_removed" --metacell_obs_col mc_truth --doublets_removed --compression none
step "run_metacells.absent" run_absent "$PYTHON" "$RUN_TOOL" \
    --cells_h5ad "$FIX/cells_nodoublets.h5ad" --gene_table "$WORKDIR/gene_table.tsv" --sample_id S1_absent \
    --prefix "$OUT/S1_absent" --metacell_obs_col mc_truth --compression none
record PASS "run_metacells" "summaries and h5ads written for four samples"

SUMMARY="$OUT/S1_starsolo_mc_summary.json"

# Pooled over the group's UMIs, on the full matrix -- not after excluding mito genes
assert "summary_qc.mito" \
    "((lambda s: (all(close(mc(s, g)['mito_pct'], F.expected_mito_pct(g)) for g in F.PROFILES),
                  {g: mc(s, g)['mito_pct'] for g in F.PROFILES}))(load(A[0])))" \
    "pooled mito % per group equals the planted profile" \
    "$SUMMARY"

assert "summary_qc.doublets" \
    "((lambda s: (all(close(mc(s, g)['doublet_pct'], F.expected_doublet_pct(g)) for g in F.PROFILES)
                  and s['sample']['doublet_mode'] == 'annotated',
                  {g: mc(s, g)['doublet_pct'] for g in F.PROFILES}))(load(A[0])))" \
    "doublet % per group equals the planted calls" \
    "$SUMMARY"

# alpha_hat comes from the shuffled raw object: a lift by position would give the
# empty droplets' 0.99
assert "summary_qc.ambient" \
    "((lambda s: (all(close(mc(s, g)['alpha_hat_mean'], F.ALPHA[g], 1e-9) for g in F.PROFILES)
                  and s['sample']['ambient_available'],
                  {g: mc(s, g)['alpha_hat_mean'] for g in F.PROFILES}))(load(A[0])))" \
    "alpha_hat lifted by barcode from the raw object" \
    "$SUMMARY"

assert "summary_qc.groups" \
    "((lambda s: ([r['id'] for r in s['metacells']] == F.GROUPS + ['__outliers__', '__excluded__']
                  and [r['pseudo'] for r in s['metacells']] == [False] * 3 + [True] * 2
                  and s['sample']['n_metacells'] == 3
                  and [r['n_cells'] for r in s['metacells']] == [F.group_sizes()[g] for g in F.PROFILES],
                  [r['id'] for r in s['metacells']]))(load(A[0])))" \
    "real metacells in order, pseudo-groups last and flagged" \
    "$SUMMARY"

assert "summary_qc.markers" \
    "((lambda s: (mc(s, 'M0.1')['markers'][0] == 'Actb' and not any(m.startswith('mt-') for r in s['metacells'] for m in r['markers']),
                  {r['id']: r['markers'] for r in s['metacells'] if not r['pseudo']}))(load(A[0])))" \
    "the most enriched gene leads; mitochondrial genes are never markers" \
    "$SUMMARY"

assert "summary_modes" \
    "((lambda r, a: (r['sample']['doublet_mode'] == 'removed' and a['sample']['doublet_mode'] == 'absent'
                     and all(x['doublet_pct'] is None for x in r['metacells'] + a['metacells'])
                     and not r['sample']['ambient_available'],
                     f\"{r['sample']['doublet_mode']} / {a['sample']['doublet_mode']}\"))(load(A[0]), load(A[1])))" \
    "removed and absent calls give null, never 0 %" \
    "$OUT/S1_removed_mc_summary.json" "$OUT/S1_absent_mc_summary.json"

assert "cells_h5ad" \
    "((lambda a: (list(a.obs['metacell'].astype(str)) == F.groups()
                  and np.allclose(a.obs['alpha_hat'].to_numpy(dtype=float), F.alpha_hat())
                  and np.array_equal(dense(a.layers['cellsweep']), F.denoised(F.counts()))
                  and np.array_equal(dense(a.X), F.counts())
                  and np.allclose(a.obs['pct_mito'], F.expected_cell_mito_pct())
                  and list(a.var['is_mito']) == [g in F.MITO_GENES for g in F.GENE_IDS],
                  'obs, layer and var follow the barcode and gene'))(read_h5ad(A[0])))" \
    "X raw, CellSweep's layer and QC on their own cells" \
    "$OUT/S1_starsolo_mc2_cells.h5ad"

assert "fingerprint" \
    "((lambda a, s, s2: (a.uns['mc_fingerprint'] == s['sample']['fingerprint'] == s2['sample']['fingerprint']
                         == qcs.fingerprint(a.obs_names[::-1], a.obs['metacell'].astype(str)[::-1], a.var_names)
                         and qcs.fingerprint(F.barcodes(), ['M1.2'] + F.groups()[1:], F.GENE_IDS) != s['sample']['fingerprint'],
                         s['sample']['fingerprint'][:23]))(read_h5ad(A[0]), load(A[1]), load(A[2])))" \
    "order-independent, stored with the cells, changed by one reassigned cell" \
    "$OUT/S1_starsolo_mc2_cells.h5ad" "$SUMMARY" "$OUT/S1_alevinfry_mc_summary.json"

# --------------------------------------------------------------------------
# The report
# --------------------------------------------------------------------------

REPORT="$WORKDIR/filtering_report.html"
step "report.run" report "$PYTHON" "$REPORT_TOOL" \
    --template "$TEMPLATE" \
    --summaries "$OUT/S1_starsolo_mc_summary.json" "$OUT/S1_alevinfry_mc_summary.json" \
    --mt_rrna_metrics "$FIX/S1_starsolo_mt_rrna_metrics.txt" \
    --antisense_metrics "$FIX/S1_starsolo_antisense_metrics.txt" \
    --max_mito_pct 20 --max_doublet_pct 50 --min_gene_umis 2 \
    --pfam_presets '^Ribosomal;^PF00125' --version test \
    --logo "$PROJECT_ROOT/assets/bca_logo.svg" --output "$REPORT"

assert "report_escape" \
    "((lambda html: (lambda blocks: (
        len(blocks) == 1 and all(json.loads(b) is not None for b in blocks)
        and '<script>alert' not in html and html.count('</script>') == 3
        and '__REPORT_DATA_PLACEHOLDER__' in json.loads(blocks[0])['universes'][0]['names']
        and 'class=\"logo\"' in html and 'src=\"data:image/svg+xml;base64,' in html
        and '__LOGO_PLACEHOLDER__' not in html,
        f'{len(blocks)} JSON block, logo inlined, {len(html)} bytes'))(
            __import__('re').findall(r'<script id=\"report-data\" type=\"application/json\">(.*?)</script>', html, __import__('re').S)))(
        open(A[0], encoding='utf-8').read()))" \
    "the payload parses, the logo is inlined; neither the script name nor the placeholder name broke out" \
    "$REPORT"

PAYLOAD_EXPR="json.loads(__import__('re').search(r'<script id=\"report-data\" type=\"application/json\">(.*?)</script>', open(A[0], encoding='utf-8').read(), __import__('re').S).group(1))"

assert "report_payload.universe" \
    "((lambda p: (len(p['universes']) == 1 and [s['universe'] for s in p['samples']] == [0, 0]
                  and F.ZERO_GENE not in p['universes'][0]['ids'] and p['universes'][0]['n_unexpressed'] == 1
                  and sorted(p['domains']) == ['PF00125', 'Ribosomal_L3'],
                  f\"{len(p['universes'])} universe of {len(p['universes'][0]['ids'])} genes\"))($PAYLOAD_EXPR))" \
    "both samples share one gene axis; the silent gene is left out" \
    "$REPORT"

assert "report_payload.featurecounts" \
    "((lambda p: (lambda fc: (
        close(fc['S1_starsolo']['mt_reads_pct'], 100 * F.FC_MT_FRACTION)
        and close(fc['S1_starsolo']['rrna_reads_pct'], 100 * F.FC_RRNA_FRACTION)
        and close(fc['S1_starsolo']['antisense_pct'], 100 * F.FC_ANTISENSE_FRACTION)
        and fc['S1_starsolo']['mt_reads_scope'] == 'called cells'
        and fc['S1_starsolo']['rrna_reads_scope'] == 'called cells'
        and not fc['S1_starsolo']['from_other_mapper']
        and fc['S1_alevinfry']['source_id'] == 'S1_starsolo' and fc['S1_alevinfry']['from_other_mapper'],
        fc))({s['id']: s['featurecounts'] for s in p['samples']}))($PAYLOAD_EXPR))" \
    "fractions shown as %, called cells preferred, alevin-fry labelled as borrowing STARsolo's" \
    "$REPORT"

# Without called-cell rows the whole-BAM value is shown and labelled so. A file from
# an earlier version is not read at all: its mtDNA row counted alignment records and
# its rRNA row had another denominator, so showing them under the new label would mislead
assert "report.featurecounts_scope" \
    "((lambda g: (lambda cells, bam: (
        g.pick_scoped({bam: '0.0100'}, cells, bam) == (1.0, 'whole BAM')
        and g.pick_scoped({cells: 'N/A', bam: '0.0200'}, cells, bam) == (2.0, 'whole BAM')
        and g.pick_scoped({'Percentage of mtDNA reads (of mapped reads)': '0.5',
                           'Percentage of rRNA reads (of uniquely mapped reads)': '0.5'}, cells, bam) == (None, None),
        'whole-BAM fallback labelled; legacy rows ignored'))(*g.FC_KEYS['mt_reads_pct']))(
        __import__('generate_filtering_report')))" \
    "whole BAM when no called-cell value; no legacy alias"

# The GeneExt alevin-fry run borrows from the STARsolo run on the same extended
# annotation, not from the standard one
assert "report.geneext_alevinfry" \
    "((lambda g: (lambda fc: (
        g.base_id('S1_geneext_alevinfry') == 'S1' and g.base_id('S1_subsampled_starsolo') == 'S1'
        and fc('S1_geneext_alevinfry')['source_id'] == 'S1_geneext_starsolo'
        and fc('S1_alevinfry')['source_id'] == 'S1_starsolo',
        f\"{fc('S1_geneext_alevinfry')['source_id']}, {fc('S1_alevinfry')['source_id']}\"))(
            lambda sid: g.featurecounts_for(sid, {'S1_starsolo': '', 'S1_geneext_starsolo': ''}, {})))(
        __import__('generate_filtering_report')))" \
    "_geneext_alevinfry pairs with _geneext_starsolo"

# In 'alevin_subsampled_starsolo' mode the GeneExt STARsolo run is the subsampled one
assert "report.geneext_subsampled" \
    "((lambda g: (lambda fc: (
        g.base_id('S1_geneext_subsampled_starsolo') == 'S1'
        and fc('S1_geneext_alevinfry')['source_id'] == 'S1_geneext_subsampled_starsolo',
        f\"{g.base_id('S1_geneext_subsampled_starsolo')}, {fc('S1_geneext_alevinfry')['source_id']}\"))(
            lambda sid: g.featurecounts_for(sid, {'S1_subsampled_starsolo': '', 'S1_geneext_subsampled_starsolo': ''}, {})))(
        __import__('generate_filtering_report')))" \
    "_geneext_alevinfry pairs with _geneext_subsampled_starsolo"

# The samplesheet's expected cells reach the sample they belong to; a sample whose row
# gives none carries null
assert "report.expected_cells" \
    "((lambda g: (lambda p: (
        p == {'S1_starsolo': 5000, 'a=b_alevinfry': 300}
        and [s['expected_cells'] for s in g.build_payload(
            [{'sample': {'id': 'S1_starsolo'}}, {'sample': {'id': 'S2_starsolo'}}], {}, {}, {}, {}, p)['samples']] == [5000, None],
        str(p)))(g.parse_expected_cells(['S1_starsolo=5000', 'a=b_alevinfry=300'])))(
        __import__('generate_filtering_report')))" \
    "id=n pairs parsed (id may hold '='), unlisted sample null"

# --------------------------------------------------------------------------
# The selection: from the page's own script when possible
#
# Choices, S1_starsolo: mito > 20 % excludes M1.2 (30 %); doublets > 70 % keeps
# M2.3 (66.7 %) and M0.1. S1_alevinfry: thresholds off, outliers kept, M0.1
# blacklisted. Genes: fewer than 2 UMIs (GENE011), domain PF00125 (GENE004), and
# Actb by name.
# --------------------------------------------------------------------------

CHOICES="$WORKDIR/choices.json"
cat >"$CHOICES" <<'EOF'
{"samples": {"S1_starsolo": {"maxMito": 20, "maxDoublet": 70},
             "S1_alevinfry": {"maxMito": null, "maxDoublet": null, "includeOutliers": true, "blacklist": ["M0.1"]}},
 "genes": {"minUmis": 2, "pfam": ["PF00125"], "manual": "Actb"}}
EOF

SELECTION="$WORKDIR/selection.json"
if [[ "$HAVE_NODE" -eq 1 ]]; then
    PANEL="$WORKDIR/summary_panel.json"
    if run_logged "$BCA_TEST_LOGDIR/metacell_filtering_report_js.log" "$NODE" "$HARNESS" "$REPORT" "$CHOICES" "$PANEL" \
            && sed -n '/^{/,$p' "$BCA_TEST_LOGDIR/metacell_filtering_report_js.log" >"$SELECTION" \
            && [[ -s "$SELECTION" ]]; then
        record PASS "report_js.run" "the page initialises, renders every tab and exports"
    else
        record FAIL "report_js.run" "see metacell_filtering_report_js.log"
        tail_log "$BCA_TEST_LOGDIR/metacell_filtering_report_js.log"
        finish_check
        exit 1
    fi

    assert "report_js.cells" \
        "((lambda s: (s['samples']['S1_starsolo']['keep_metacells'] == ['M0.1', 'M2.3']
                      and s['samples']['S1_alevinfry']['keep_metacells'] == ['M1.2', 'M2.3', '__outliers__']
                      and not any('whitelist' in v for v in s['samples'].values()),
                      {k: v['keep_metacells'] for k, v in s['samples'].items()}))(load(A[0])))" \
        "thresholds, the blacklist and the outlier toggle; no whitelist in the export" \
        "$SELECTION"

    assert "report_js.genes" \
        "((lambda s: (s['samples']['S1_starsolo']['excluded_genes'] == {'GENE001': ['manual'], 'GENE004': ['pfam:PF00125'], 'GENE011': ['min_umi']}
                      and s['schema_version'] == 2
                      and s['samples']['S1_starsolo']['gene_rules']['min_total_umis'] == 2
                      and s['samples']['S1_starsolo']['fingerprint'].startswith('sha256:'),
                      s['samples']['S1_starsolo']['excluded_genes']))(load(A[0])))" \
        "minimum UMIs, ticked domain and a gene given by name" \
        "$SELECTION"

    # The summary panel shows the first sample, S1_alevinfry: of 100 cells, M0.1 (30,
    # blacklisted) and __excluded__ (5) go; M1.2, M2.3 and the 5 outliers stay. Of the
    # 10 expressed genes, GENE001, GENE004 and GENE011 go.
    assert "report_js.summary" \
        "((lambda p: (all(x in p['stats'] for x in ('Called cells', '>100<', 'Metacells', 'Outlier cells'))
                      and not any(x in p['stats'] for x in ('Doublet calls', 'Ambient RNA'))
                      and all(x in p['settings'] for x in ('Cell level', 'Gene level', 'Blacklist', '>1 metacell<',
                                                            'PFAM domains', '>1 domain<', '>≥ 2<'))
                      and not any(x in p['settings'] for x in ('excludes', 'Outlier cells', 'Ungrouped cells'))
                      and all(x in p['totals'] for x in ('2 included<span class=\"sub\">1 excluded',
                                                          '65 included<span class=\"sub\">35 excluded',
                                                          '7 included<span class=\"sub\">3 excluded',
                                                          '% of UMIs kept'))
                      and 'UMIs in excluded genes' not in p['totals'],
                      'statistics, both settings tables and running totals'))(load(A[0])))" \
        "both columns of the summary panel, with totals across both levels" \
        "$PANEL"

    # A dragged line sets its threshold from x alone (rounded to 0.1 %); the rank plot's
    # line sets the minimum in UMIs; a bar click at log10(UMIs + 1) = 2 means 99 UMIs.
    assert "report_js.drag" \
        "((lambda d: (d['mito'] == 12.4 and d['rank'] == 50 and d['click'] == 99, d))(load(A[0])['drag']))" \
        "dragging the threshold lines and clicking a bar move the thresholds" \
        "$PANEL"

    # Every sample has its own settings at both levels. With S1_starsolo's gene level pending,
    # the selection holds S1_alevinfry alone and the export table shows its genes kept as
    # pending. Opening the gene tab with S1_starsolo selected applies its gene level only.
    # Applying S1_alevinfry's settings to all copies the mito threshold (33) and the minimum
    # UMIs (7) and applies both levels, while S1_starsolo keeps its own blacklists. A
    # version-1 selection (one set of gene rules for the run) reads back the same choices.
    assert "report_js.flow" \
        "((lambda f: (f['with_pending'] == ['S1_alevinfry'] and f['export_pending'] == 1
                      and f['visited'] == {'cells': False, 'genes': True}
                      and f['applied'] == {'maxMito': 33, 'minUmis': 7, 'set': {'cells': True, 'genes': True}, 'kept_blacklist': True}
                      and f['v1_problems'] == [] and f['v1_same'],
                      f))(load(A[0])['flow']))" \
        "per-sample levels: a pending level leaves its sample out, a visit applies it, Apply to all keeps blacklists, version 1 reads back" \
        "$PANEL"
else
    record SKIP "report_js" "no JavaScript runtime (pass --node, or set BCA_NODE); selection built in Python"
    "$PYTHON" - "$SUMMARY" "$OUT/S1_alevinfry_mc_summary.json" "$SELECTION" <<'PYEOF'
import json, sys
s1, s2, out = sys.argv[1:4]
fp = lambda p: json.load(open(p))["sample"]["fingerprint"]
excluded = {"GENE001": ["manual"], "GENE004": ["pfam:PF00125"], "GENE011": ["min_umi"]}
rules = {"min_total_umis": 2}
json.dump({"schema": "bca_metacell_selection", "schema_version": 2,
           "samples": {"S1_starsolo": {"fingerprint": fp(s1), "keep_metacells": ["M0.1", "M2.3"], "excluded_genes": excluded, "gene_rules": rules},
                       "S1_alevinfry": {"fingerprint": fp(s2), "keep_metacells": ["M1.2", "M2.3", "__outliers__"], "excluded_genes": excluded, "gene_rules": rules}}},
          open(out, "w"))
PYEOF
fi

# Version 2 carries each sample's gene rules; a version-1 selection's one set still applies
assert "apply_gene_rules" \
    "((lambda a: (a.gene_rules_for({'gene_rules': {'min_total_umis': 1}}, {'gene_rules': {'min_total_umis': 5}}) == {'min_total_umis': 5}
                  and a.gene_rules_for({'gene_rules': {'min_total_umis': 1}}, {}) == {'min_total_umis': 1},
                  'the sample own rules first, else the run rules'))(__import__('apply_metacell_filter')))" \
    "per-sample gene rules (version 2), the run's (version 1) otherwise"

# --------------------------------------------------------------------------
# The apply step
# --------------------------------------------------------------------------

FINAL="$WORKDIR/final"
step "apply.run" apply "$PYTHON" "$APPLY_TOOL" \
    --cells_h5ad "$OUT/S1_starsolo_mc2_cells.h5ad" --selection "$SELECTION" \
    --sample_id S1_starsolo --outdir "$FINAL" --compression none
FLAGGED="$WORKDIR/final_flag"
step "apply.run_flag" apply_flag "$PYTHON" "$APPLY_TOOL" \
    --cells_h5ad "$OUT/S1_starsolo_mc2_cells.h5ad" --selection "$SELECTION" \
    --sample_id S1_starsolo --outdir "$FLAGGED" --gene_mode flag --compression none

KEPT_EXPR="[i for i, g in enumerate(F.groups()) if g in ('M0.1', 'M2.3')]"
DROPPED_EXPR="['GENE001', 'GENE004', 'GENE011', F.ZERO_GENE]"

assert "apply_keep" \
    "((lambda m, keep: (lambda gi: (
        m[1] == [F.barcodes()[i] for i in keep]
        and np.array_equal(dense(m[0]), F.counts()[np.ix_(keep, gi)].T),
        f'{len(m[1])} cells x {m[0].shape[0]} genes'))(
            [j for j, g in enumerate(F.GENE_IDS) if g not in $DROPPED_EXPR]))(
        __import__('mtx_io').read_triplet(A[0] + '/matrix.mtx', A[0] + '/barcodes.tsv', A[0] + '/features.tsv'), $KEPT_EXPR))" \
    "the kept metacells' cells in order, every count unchanged" \
    "$FINAL"

assert "apply_genes.drop" \
    "((lambda a, s: (sorted(set(F.GENE_IDS) - set(a.var_names)) == sorted($DROPPED_EXPR)
                     and s['genes_excluded'] == 4 and s['cells_kept'] == 60 and s['groups_kept'] == 2,
                     f\"dropped {sorted(set(F.GENE_IDS) - set(a.var_names))}\"))(
        read_h5ad(A[0] + '/S1_starsolo_final.h5ad'), load(A[0] + '/S1_starsolo_filter_summary.json')))" \
    "listed genes and the unexpressed one (by the min-UMI rule) are dropped" \
    "$FINAL"

assert "apply_genes.flag" \
    "((lambda a: (list(a.var_names) == F.GENE_IDS
                  and sorted(a.var_names[~a.var['pass_filter'].to_numpy(dtype=bool)]) == sorted($DROPPED_EXPR)
                  and a.var.loc[F.ZERO_GENE, 'exclusion_reasons'] == 'min_umi',
                  'full gene axis, decisions in var'))(read_h5ad(A[0] + '/S1_starsolo_final.h5ad')))" \
    "flag mode keeps every gene and records why" \
    "$FLAGGED"

assert "apply_metacell_sums" \
    "((lambda m: (list(m.obs_names) == ['M0.1', 'M2.3']
                  and np.allclose(dense(m.X)[0].sum(), F.counts()[[i for i, g in enumerate(F.groups()) if g == 'M0.1']][:, [j for j, g in enumerate(F.GENE_IDS) if g not in $DROPPED_EXPR]].sum()),
                  list(m.obs_names)))(read_h5ad(A[0] + '/S1_starsolo_final_metacells.h5ad')))" \
    "metacell x gene sums of the kept cells" \
    "$FINAL"

# A guard that passes silently is worse than none, so each must exit with its code
"$PYTHON" - "$SELECTION" "$WORKDIR" <<'PYEOF'
import json, sys
sel = json.load(open(sys.argv[1]))
stale = json.loads(json.dumps(sel)); stale["samples"]["S1_starsolo"]["fingerprint"] = "sha256:" + "0" * 64
json.dump(stale, open(sys.argv[2] + "/sel_stale.json", "w"))
unknown = json.loads(json.dumps(sel)); unknown["samples"]["S1_starsolo"]["keep_metacells"].append("M99.9")
json.dump(unknown, open(sys.argv[2] + "/sel_unknown.json", "w"))
other = json.loads(json.dumps(sel)); other["schema"] = "something_else"
json.dump(other, open(sys.argv[2] + "/sel_schema.json", "w"))
PYEOF

GUARDS_OK=1
for guard in stale:S1_starsolo unknown:S1_starsolo schema:S1_starsolo stale:S1_missing; do
    kind="${guard%%:*}" sid="${guard#*:}"
    set +e
    run_logged "$BCA_TEST_LOGDIR/metacell_filtering_guard_${kind}_${sid}.log" "$PYTHON" "$APPLY_TOOL" \
        --cells_h5ad "$OUT/S1_starsolo_mc2_cells.h5ad" --selection "$WORKDIR/sel_${kind}.json" \
        --sample_id "$sid" --outdir "$WORKDIR/guard_${kind}_${sid}" --compression none
    rc=$?
    set -e
    if [[ "$rc" -ne 2 ]]; then
        record FAIL "apply_guards.$kind.$sid" "exit $rc, expected 2"
        GUARDS_OK=0
    elif [[ -e "$WORKDIR/guard_${kind}_${sid}/matrix.mtx" ]]; then
        record FAIL "apply_guards.$kind.$sid" "wrote a matrix despite refusing"
        GUARDS_OK=0
    fi
done
[[ "$GUARDS_OK" -eq 1 ]] && record PASS "apply_guards" "stale fingerprint, unknown metacell, foreign schema, absent sample: exit 2, nothing written"

# --------------------------------------------------------------------------
# Case: mc2_smoke -- the real grouping, where metacells is installed
# --------------------------------------------------------------------------

if ! "$PYTHON" -c 'import metacells' >/dev/null 2>&1; then
    record SKIP "mc2_smoke" "metacells is not importable (it runs in the METACELL2 environment)"
else
    SMOKE_OK=1
    for run in a b; do
        if ! run_logged "$BCA_TEST_LOGDIR/metacell_filtering_mc2_${run}.log" "$PYTHON" "$RUN_TOOL" \
                --cells_h5ad "$FIX/mc2_cells.h5ad" --gene_table "$FIX/mc2_gene_table.tsv" --sample_id mc2 \
                --prefix "$WORKDIR/mc2_${run}" --target_metacell_size 30 --min_cell_umis 50 --min_cells 100 \
                --random_seed 1 --cpus 2 --compression none; then
            record FAIL "mc2_smoke.run_$run" "see metacell_filtering_mc2_${run}.log"
            SMOKE_OK=0
        fi
    done
    if [[ "$SMOKE_OK" -eq 1 ]]; then
        # Metacells of a planted type should not straddle types: each holds one type's cells
        assert "mc2_smoke" \
            "((lambda a, b, c: (a['sample']['status'] == 'ok' and a['sample']['n_metacells'] >= F.SMOKE_TYPES
                                and a['sample']['fingerprint'] == b['sample']['fingerprint']
                                and all(len(set(g)) == 1 for m, g in
                                        pd.Series(c.obs['truth'].astype(str).to_numpy()).groupby(c.obs['metacell'].astype(str).to_numpy())
                                        if not m.startswith('__'))
                                and all(x['x'] is not None for x in a['metacells'] if not x['pseudo']),
                                f\"{a['sample']['n_metacells']} pure metacells with a UMAP, same fingerprint twice\"))(
                load(A[0]), load(A[1]), read_h5ad(A[2])))" \
            "Metacell2 recovers the planted types, the UMAP is computed, the same seed gives the same assignment" \
            "$WORKDIR/mc2_a_mc_summary.json" "$WORKDIR/mc2_b_mc_summary.json" "$WORKDIR/mc2_a_mc2_cells.h5ad"
    fi

    if run_logged "$BCA_TEST_LOGDIR/metacell_filtering_mc2_tiny.log" "$PYTHON" "$RUN_TOOL" \
            --cells_h5ad "$FIX/cells.h5ad" --gene_table "$WORKDIR/gene_table.tsv" --sample_id tiny \
            --prefix "$WORKDIR/tiny" --min_cells 100000 --compression none; then
        assert "mc2_smoke.skipped" \
            "((lambda s: (s['sample']['status'] == 'skipped' and s['sample']['reason'] and not __import__('os').path.exists(A[1]),
                          s['sample']['reason']))(load(A[0])))" \
            "too few cells: a 'skipped' summary, exit 0, no h5ad" \
            "$WORKDIR/tiny_mc_summary.json" "$WORKDIR/tiny_mc2_cells.h5ad"
    else
        record FAIL "mc2_smoke.skipped" "a sample with too few cells must exit 0"
    fi
fi

finish_check
