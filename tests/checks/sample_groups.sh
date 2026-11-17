#!/usr/bin/env bash
# description: Check the samplesheet's qc_report and filtering_report groups: pooling, report and split
#
# A filtering_report group pools its samples' cell-called objects before Metacell2, and
# everything after reads the pool: the summary, the report's reads and read metrics,
# and the filtered matrices, which are split back per sample. Each link can go wrong
# without raising: a cell given the other sample's CellSweep values, a read count taken
# from the wrong CellReads.stats, a barcode that keeps its '<sample>:' prefix. The
# qc_report group only pools numbers in the dashboard, where a rate averaged instead of
# recomputed from the sums reads plausibly.
#
# tests/lib/make_group_fixture.py builds two samples on make_metacell_fixture.py's,
# with every value of the second shifted so its origin can be told apart.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

MERGE_TOOL="$PROJECT_ROOT/bin/merge_group_h5ad.py"
GENE_TABLE_TOOL="$PROJECT_ROOT/bin/build_gene_table.py"
RUN_TOOL="$PROJECT_ROOT/bin/run_metacells.py"
REPORT_TOOL="$PROJECT_ROOT/bin/generate_filtering_report.py"
APPLY_TOOL="$PROJECT_ROOT/bin/apply_metacell_filter.py"
TEMPLATE="$PROJECT_ROOT/bin/filtering_report.html"
FIXTURE_GEN="$TESTS_DIR/lib/make_group_fixture.py"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh sample_groups [-- OPTIONS]
       tests/checks/sample_groups.sh [OPTIONS]

Validate the pooling of samplesheet groups.

Options:
  --python PATH   Python interpreter to use (default: \$BCA_PYTHON, else
                  python3, else python). Needs numpy, pandas, scipy, anndata, h5py.
  --keep          Keep the generated fixtures and outputs for inspection.
  -h, --help      Show this message.

Cases:
  merge_cells           cells named '<sample>:<barcode>', obs['sample'], members in uns
  merge_cellsweep       each sample's CellSweep values lifted from its own raw object;
                        ambient_hat averaged over the cells, each sample's kept
  merge_genes           a sample counted against other genes is refused
  merge_shared          one library under two names is warned about
  summary_members       the summary lists the members, each cell's member and each
                        metacell's cells per member
  report_reads          each cell's reads from its own sample's CellReads.stats
  report_metrics        read metrics averaged over the samples' cells, each kept
  apply_members         the filtered group split per sample under the original barcodes
  dashboard_pool        sums, rates from the sums, weighted means and pooled medians
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --python)  PYTHON="${2:?--python needs a value}"; shift 2 ;;
        --keep)    KEEP=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

log_header "sample_groups"

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

WORKDIR="$BCA_TEST_LOGDIR/sample_groups"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# assert NAME PYTHON_EXPR DETAIL [PATH...] -- as in metacell_filtering.sh, with the group
# fixture importable as G
assert() {
    local name="$1" expr="$2" detail="${3:-}"; shift 3 || shift $#
    local out
    if out="$("$PYTHON" -c "
import sys, json, warnings
warnings.simplefilter('ignore')
LIB, BIN, W, A = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
sys.path.insert(0, LIB)
sys.path.insert(0, BIN)
import numpy as np, pandas as pd
import metacell_qc_summary as qcs
import make_metacell_fixture as F
import make_group_fixture as G

def load(path):
    with open(path, encoding='utf-8') as fh:
        return json.load(fh)

def read_h5ad(path):
    import anndata as ad
    return ad.read_h5ad(path)

def dense(m):
    return np.asarray(m.todense()) if hasattr(m, 'todense') else np.asarray(m)

def payload(path):
    import re
    html = open(path, encoding='utf-8').read()
    return json.loads(re.search(r'<script id=\"report-data\" type=\"application/json\">(.*?)</script>', html, re.S).group(1))

ok, msg = ($expr)
print(('PASS' if ok else 'FAIL') + '\t' + str(msg))
" "$TESTS_DIR/lib" "$PROJECT_ROOT/bin" "$WORKDIR" "$@" 2>&1 | tail -1)"; then
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

step() {
    local name="$1" logname="$2"; shift 2
    if run_logged "$BCA_TEST_LOGDIR/sample_groups_${logname}.log" "$@"; then
        return 0
    fi
    record FAIL "$name" "see sample_groups_${logname}.log"
    tail_log "$BCA_TEST_LOGDIR/sample_groups_${logname}.log"
    finish_check
    exit 1
}

# --------------------------------------------------------------------------
# Fixture, pooled
# --------------------------------------------------------------------------

FIX="$WORKDIR/fixture"
step "fixture" fixture "$PYTHON" "$FIXTURE_GEN" "$FIX"
step "merge.help" merge_help "$PYTHON" "$MERGE_TOOL" --help

POOLED="$WORKDIR/G_starsolo_group_cells.h5ad"
step "merge.run" merge "$PYTHON" "$MERGE_TOOL" --group_id G_starsolo \
    --member S1_starsolo S1 "$FIX/cells_S1.h5ad" "$FIX/raw.h5ad" \
    --member S2_starsolo S2 "$FIX/cells_S2.h5ad" "$FIX/raw_S2.h5ad" \
    --output "$POOLED" --members_tsv "$WORKDIR/G_starsolo_group_members.tsv" --compression none
record PASS "merge" "two samples pooled"

assert "merge_cells" \
    "((lambda a: (list(a.obs_names) == ['S1:' + b for b in F.barcodes()] + ['S2:' + b for b in G.s2_barcodes()]
                  and list(a.obs['sample'].astype(str)) == ['S1'] * len(F.barcodes()) + ['S2'] * len(F.barcodes())
                  and list(a.obs['sample_id'].astype(str)) == ['S1_starsolo'] * len(F.barcodes()) + ['S2_starsolo'] * len(F.barcodes())
                  and list(a.uns['group_members']['sample']) == ['S1', 'S2']
                  and np.array_equal(dense(a.X), np.vstack([F.counts(), F.counts()]))
                  and list(a.obs['doublet_status'].astype(str)) == F.doublet_status() * 2,
                  f'{a.n_obs} cells, {a.obs.columns.size} obs columns'))(read_h5ad(A[0])))" \
    "cells prefixed with their sample, counts and doublet calls carried over" \
    "$POOLED"

assert "merge_cellsweep" \
    "((lambda a: (np.allclose(a.obs['alpha_hat'].to_numpy(dtype=float),
                              F.alpha_hat() + [x + G.ALPHA_SHIFT for x in F.alpha_hat()])
                  and np.array_equal(dense(a.layers['cellsweep']), np.vstack([F.denoised(F.counts())] * 2))
                  and np.allclose(a.var['ambient_hat'], [F.ambient_hat(j) * (1 + G.AMBIENT_SCALE) / 2 for j in range(len(F.GENE_IDS))])
                  and np.allclose(a.var['ambient_hat:S2'], [F.ambient_hat(j) * G.AMBIENT_SCALE for j in range(len(F.GENE_IDS))]),
                  'alpha_hat, the denoised layer and ambient_hat per sample'))(read_h5ad(A[0])))" \
    "each sample's CellSweep values from its own raw object; ambient_hat averaged over cells" \
    "$POOLED"

if run_logged "$BCA_TEST_LOGDIR/sample_groups_merge_genes.log" "$PYTHON" "$MERGE_TOOL" --group_id BAD \
        --member S1_starsolo S1 "$FIX/cells_S1.h5ad" --member S2_starsolo S2 "$FIX/cells_bad.h5ad" \
        --output "$WORKDIR/bad.h5ad" --compression none; then
    record FAIL "merge_genes" "a sample counted against other genes was pooled"
elif grep -q "share one reference annotation" "$BCA_TEST_LOGDIR/sample_groups_merge_genes.log"; then
    record PASS "merge_genes" "refused, naming the genes that differ"
else
    record FAIL "merge_genes" "failed without saying why; see sample_groups_merge_genes.log"
fi

step "merge.shared" merge_shared "$PYTHON" "$MERGE_TOOL" --group_id TWICE \
    --member S1_starsolo S1 "$FIX/cells_S1.h5ad" --member S1b_starsolo S1b "$FIX/cells_S1.h5ad" \
    --output "$WORKDIR/twice.h5ad" --compression none
if grep -q "one library sequenced twice" "$BCA_TEST_LOGDIR/sample_groups_merge_shared.log" \
   && ! grep -q "one library sequenced twice" "$BCA_TEST_LOGDIR/sample_groups_merge.log"; then
    record PASS "merge_shared" "warned for the same cells twice, not for $(printf '%s' 2) shared barcodes"
else
    record FAIL "merge_shared" "the shared-barcode warning is missing, or raised for unrelated samples"
fi

# --------------------------------------------------------------------------
# Summary, report and the split
# --------------------------------------------------------------------------

step "gene_table" gene_table "$PYTHON" "$GENE_TABLE_TOOL" \
    --gtf "$FIX/genes.gtf" --mt_contig chrM MT --rrna_pattern rRNA \
    --output "$WORKDIR/gene_table.tsv" --stats "$WORKDIR/gene_table_stats.json"
step "run_metacells" run "$PYTHON" "$RUN_TOOL" \
    --cells_h5ad "$POOLED" --gene_table "$WORKDIR/gene_table.tsv" --sample_id G_starsolo \
    --mapping_method fixture --prefix "$WORKDIR/G_starsolo" --metacell_obs_col mc_truth --compression none
SUMMARY="$WORKDIR/G_starsolo_mc_summary.json"

assert "summary_members" \
    "((lambda s: (lambda n: (
        [(m['sample'], m['id'], m['n_cells']) for m in s['sample']['members']] == [('S1', 'S1_starsolo', n), ('S2', 'S2_starsolo', n)]
        and s['cells']['member'] == [0] * n + [1] * n
        and all(r['member_counts'] == [F.group_sizes()[r['id']]] * 2 for r in s['metacells'])
        and s['sample']['n_cells'] == 2 * n,
        f\"{len(s['metacells'])} metacells, members {[m['sample'] for m in s['sample']['members']]}\"))(len(F.barcodes())))(load(A[0])))" \
    "members listed, each cell's member, each metacell's cells per member" \
    "$SUMMARY"

REPORT="$WORKDIR/filtering_report.html"
step "report" report "$PYTHON" "$REPORT_TOOL" --template "$TEMPLATE" --summaries "$SUMMARY" \
    --mt_rrna_metrics "$FIX/S1_starsolo_mt_rrna_metrics.txt" "$FIX/S2_starsolo_mt_rrna_metrics.txt" \
    --antisense_metrics "$FIX/S1_starsolo_antisense_metrics.txt" "$FIX/S2_starsolo_antisense_metrics.txt" \
    --cellreads_stats "S1_starsolo=$FIX/S1_starsolo_CellReads.stats" "S2_starsolo=$FIX/S2_starsolo_CellReads.stats" \
    --expected_cells G_starsolo=500 --version test --output "$REPORT"

assert "report_reads" \
    "((lambda p: (lambda c: (c['reads'] == F.cell_reads() + G.s2_reads()
                             and c['reads_source'] == 'S1_starsolo, S2_starsolo'
                             and p['samples'][0]['expected_cells'] == 500,
                             c['reads_source']))(p['samples'][0]['cells']))(payload(A[0])))" \
    "each cell's reads from its own sample, under its original barcode" \
    "$REPORT"

assert "report_metrics" \
    "((lambda fc: (abs(fc['mt_reads_pct'] - 100 * (F.FC_MT_FRACTION + G.S2_FRACTIONS['mt']) / 2) < 1e-2
                   and abs(fc['antisense_pct'] - 100 * (F.FC_ANTISENSE_FRACTION + G.S2_FRACTIONS['antisense']) / 2) < 1e-2
                   and fc['mt_reads_scope'] == 'called cells'
                   and [m['sample'] for m in fc['members']] == ['S1', 'S2'],
                   fc['source_id']))(payload(A[0])['samples'][0]['featurecounts']))" \
    "the samples' percentages averaged over their cells, each kept" \
    "$REPORT"

SELECTION="$WORKDIR/selection.json"
"$PYTHON" - "$SUMMARY" "$SELECTION" <<'PY'
import json, sys
s = json.load(open(sys.argv[1]))
json.dump({"schema": "bca_metacell_selection", "schema_version": 1, "samples": {"G_starsolo": {
    "fingerprint": s["sample"]["fingerprint"], "keep_metacells": ["M0.1"],
    "excluded_cells": {}, "excluded_genes": {}}}}, open(sys.argv[2], "w"))
PY
FINAL="$WORKDIR/G_starsolo_final"
step "apply" apply "$PYTHON" "$APPLY_TOOL" --cells_h5ad "$WORKDIR/G_starsolo_mc2_cells.h5ad" \
    --selection "$SELECTION" --sample_id G_starsolo --outdir "$FINAL" --compression none

assert "apply_members" \
    "((lambda keep: (lambda a1, a2, b1, summ: (
        list(a1.obs_names) == [F.barcodes()[i] for i in keep]
        and list(a2.obs_names) == [G.s2_barcodes()[i] for i in keep]
        and b1 == [F.barcodes()[i] for i in keep]
        and np.array_equal(dense(a2.X), F.counts()[keep])
        and summ['members_cells_kept'] == {'S1': len(keep), 'S2': len(keep)},
        f\"{len(keep)} cells per sample\"))(
            read_h5ad(A[0] + '/members/S1/S1_final.h5ad'), read_h5ad(A[0] + '/members/S2/S2_final.h5ad'),
            open(A[0] + '/members/S1/barcodes.tsv').read().split(),
            load(A[0] + '/G_starsolo_filter_summary.json')))(
        [i for i, g in enumerate(F.groups()) if g == 'M0.1']))" \
    "each sample's kept cells on their own, under the original barcodes" \
    "$FINAL"

# --------------------------------------------------------------------------
# qc_report: the dashboard's pooled row
# --------------------------------------------------------------------------

assert "dashboard_pool" \
    "((lambda D: (lambda g: (
        g['pooled']['n_reads'] == '3,000' and g['pooled']['pct_unique'] == '80.00%'
        and g['pooled']['n_cells'] == '300' and g['pooled']['mean_reads'] == '10'
        and g['pooled']['noise'] == '15.00%' and g['pooled']['median_umis'] == '30'
        and g['pooled']['reads_07'] == '3.0 M' and 'saturation' in g['approx']
        and 'median_umis' not in g['approx'] and g['members'] == ['A_starsolo', 'B_starsolo'],
        g['pooled']))(D.build_group_payload([
            dict(id='A_starsolo', base='A', suffix='_starsolo', source='starsolo', n_input=1000.0, n_unique=900.0,
                 n_cells=100.0, saturation=50.0, noise=10.0, intronic=None, rrna_mapped=None, rrna_unique=None,
                 mtdna_mapped=None, mtdna_cells=None, reads_07=1e6, median_umis=10.0, median_genes=5.0,
                 cell_umis=[10] * 100),
            dict(id='B_starsolo', base='B', suffix='_starsolo', source='starsolo', n_input=2000.0, n_unique=1500.0,
                 n_cells=200.0, saturation=70.0, noise=18.0, intronic=None, rrna_mapped=None, rrna_unique=None,
                 mtdna_mapped=None, mtdna_cells=None, reads_07=2e6, median_umis=30.0, median_genes=9.0,
                 cell_umis=[30] * 200),
        ], {'A': {'qc_report': 'grp'}, 'B': {'qc_report': 'grp'}})[0]))(__import__('generate_dashboard')))" \
    "reads and cells summed, rates from the sums, noise weighted by unique reads, median over all cells"

finish_check
