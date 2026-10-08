#!/usr/bin/env bash
# description: Check that the mtDNA / rRNA metrics count reads, not alignments, per library and per cell
#
# A BAM holds one record per alignment, so a count of records over-weights every
# multimapper by its NH. Nothing about the result looks wrong: the percentages stay
# plausible, they are just of the wrong thing. The fixture below is a hand-written
# SAM in which every number can be derived by hand -- a nuclear multimapper with a
# secondary alignment on chrM, an mtDNA multimapper on mitochondrial rRNA, reads
# sense, antisense and outside the mitochondrial genes, an rRNA read on an added
# rRNA contig, and reads from a called and an uncalled barcode -- and the cases
# assert those exact numbers for calculate_read_metrics.py (CALC_READ_METRICS) and
# per-cell_images.py (PERCELL_METRICS). The antisense cases read a hand-written
# STARsolo CellReads.stats instead of the BAM.
#
# Reads (all 50M):
#   r1   chrM   150  +  NH1  AAAA  sense to MT-CO1
#   r2   chrM   160  -  NH1  AAAA  antisense to MT-CO1
#   r3   chrM   700  +  NH1  CCCC  outside every mitochondrial gene
#   r4   chrM  1350  +  NH2  AAAA  on MT-RNR2 (Mt_rRNA); secondary on chr1
#   r5   chr1  3100  +  NH3  AAAA  protein coding; secondaries on chrM and chr1
#   r6   chr1  1100  +  NH1  AAAA  nuclear rRNA gene
#   r7   rDNA   100  +  NH1  CCCC  added rRNA contig
#   r8   unmapped            AAAA
#   r9   chr1  3200  +  NH1  GGGG
#   r10  chr1  3300  +  NH1  AAAA
#
# Needs samtools (to build the fixture) and a Python with pysam; the per-cell case
# also needs pandas, scipy and matplotlib in that Python.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

METRICS_TOOL="$PROJECT_ROOT/bin/calculate_read_metrics.py"
PERCELL_TOOL="$PROJECT_ROOT/bin/per-cell_images.py"

PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh mt_rrna_metrics [-- OPTIONS]
       tests/checks/mt_rrna_metrics.sh [OPTIONS]

Validate the read-level mtDNA / rRNA metrics against a hand-counted BAM.

Options:
  --python PATH   Python interpreter (default: \$BCA_PYTHON, else python3). Needs
                  pysam; the percell case also pandas, scipy and matplotlib.
  --keep          Keep the generated fixture and outputs for inspection.
  -h, --help      Show this message.

Cases:
  library         mapped, unique and multimapped reads; rRNA and mtDNA of each
  origin          mtDNA reads sense, antisense and outside mitochondrial genes
  unstranded      without a strand only "outside" is reported
  called_cells    the same metrics restricted to the filtered barcodes
  no_cells        without a barcode list the called-cell rows are N/A
  antisense       antisense share of gene reads from CellReads.stats, library and cells
  antisense_unstranded  no antisense file for an unstranded run
  barcode_reads   per-barcode reads; the called barcodes sum to the called-cell rows
  percell         per-cell mapped reads, mtDNA %, rRNA %, intronic % and unspliced %,
                  from the per-barcode reads of the stranded run
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

log_header "mt_rrna_metrics"

if [[ -z "$PYTHON" ]]; then
    if have_cmd python3; then PYTHON="python3"; elif have_cmd python; then PYTHON="python"; fi
fi
if ! have_cmd samtools || [[ -z "$PYTHON" ]] || ! "$PYTHON" -c 'import pysam' >/dev/null 2>&1; then
    record SKIP "tools" "samtools and a Python with pysam are needed"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/mt_rrna_metrics"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# --------------------------------------------------------------------------
# Fixture
# --------------------------------------------------------------------------

FIX="$WORKDIR/fixture"
mkdir -p "$FIX"

# The added rRNA reference carries no biotype at all: its contig counts as rRNA
# as a whole. The pipeline's ref_gtf arrives merged with it, so it is appended.
printf 'rDNA\tadded\texon\t1\t500\t.\t+\t.\tgene_id "rDNA"; transcript_id "rDNA.1";\n' > "$FIX/added.gtf"
{
    gtf() { printf '%s\ttest\t%s\t%s\t%s\t.\t%s\t.\tgene_id "%s"; transcript_id "%s.1"; gene_biotype "%s";\n' "$@"; }
    for ft in gene exon; do
        gtf chr1 "$ft" 1001 1500 + RNA45S RNA45S rRNA
        gtf chr1 "$ft" 3001 4000 + GENE1 GENE1 protein_coding
        gtf chrM "$ft" 101 600 + MT-CO1 MT-CO1 protein_coding
        gtf chrM "$ft" 801 1200 - MT-ND6 MT-ND6 protein_coding
        gtf chrM "$ft" 1301 1500 + MT-RNR2 MT-RNR2 Mt_rRNA
    done
    cat "$FIX/added.gtf"
} > "$FIX/ref.gtf"

{
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    printf '@SQ\tSN:chr1\tLN:10000\n'
    printf '@SQ\tSN:chrM\tLN:2000\n'
    printf '@SQ\tSN:rDNA\tLN:500\n'
    sam() {  # qname flag rname pos nh cb
        if [[ "$2" == 4 ]]; then
            printf '%s\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tNH:i:0\tCB:Z:%s\n' "$1" "$6"
        else
            printf '%s\t%s\t%s\t%s\t255\t50M\t*\t0\t0\t*\t*\tNH:i:%s\tCB:Z:%s\n' "$@"
        fi
    }
    sam r1  0   chrM 150  1 AAAA
    sam r2  16  chrM 160  1 AAAA
    sam r3  0   chrM 700  1 CCCC
    sam r4  0   chrM 1350 2 AAAA
    sam r4  256 chr1 5000 2 AAAA
    sam r5  0   chr1 3100 3 AAAA
    sam r5  256 chrM 200  3 AAAA
    sam r5  256 chr1 6000 3 AAAA
    sam r6  0   chr1 1100 1 AAAA
    sam r7  0   rDNA 100  1 CCCC
    sam r8  4   '*'  0    0 AAAA
    sam r9  0   chr1 3200 1 GGGG
    sam r10 0   chr1 3300 1 AAAA
} > "$FIX/reads.sam"

samtools sort -o "$FIX/sample.bam" "$FIX/reads.sam" 2>/dev/null
samtools index "$FIX/sample.bam"

# Called cells: AAAA, plus one barcode with no reads at all
printf 'AAAA\nTTTT\n' > "$FIX/barcodes.tsv"

# STARsolo read flags for the antisense metrics. Library (every row): exonic 68,
# intronic 23, exonicAS 12, intronicAS 6. Called cells (AAAA only): 60, 20, 10, 5.
{
    printf 'CB\tgenomeU\tgenomeM\texonic\tintronic\texonicAS\tintronicAS\tmito\n'
    printf 'CBnotInPasslist\t10\t0\t3\t2\t1\t0\t0\n'
    printf 'AAAA\t100\t5\t60\t20\t10\t5\t3\n'
    printf 'CCCC\t10\t0\t5\t1\t1\t1\t0\n'
    printf 'GGGG\t0\t0\t0\t0\t0\t0\t0\n'
} > "$FIX/CellReads.stats"

# run_metrics NAME OUTDIR ARGS... -- runs the tool in its own directory, so each
# case keeps its outputs apart
run_metrics() {
    local name="$1" out="$2"; shift 2
    mkdir -p "$out"
    if ! (cd "$out" && run_logged "$BCA_TEST_LOGDIR/mt_rrna_metrics_${name}.log" \
            "$PYTHON" "$METRICS_TOOL" "$@" --bam "$FIX/sample.bam" --gtf "$FIX/ref.gtf" \
                --mt-contig "chrM M MT" --cellreads "$FIX/CellReads.stats" \
                --out-mt-rrna metrics.csv --out-antisense antisense.csv \
                --out-barcode-reads barcode_reads.tsv.gz); then
        record FAIL "$name" "calculate_read_metrics.py failed, see mt_rrna_metrics_${name}.log"
        tail_log "$BCA_TEST_LOGDIR/mt_rrna_metrics_${name}.log"
        return 1
    fi
}

# value FILE METRIC -- the value of one metric row; names are quoted in the file
value() {
    awk -F'",' -v k="\"$2" '$1 == k { v = $2; gsub(/"/, "", v); print v; found = 1 }
                            END { if (!found) print "<missing>" }' "$1"
}

# expect CASE FILE "METRIC=VALUE"... -- PASS when every metric matches exactly
expect() {
    local name="$1" file="$2"; shift 2
    local bad="" pair metric want got
    for pair in "$@"; do
        metric="${pair%=*}"; want="${pair##*=}"
        got="$(value "$file" "$metric")"
        [[ "$got" == "$want" ]] || bad="${bad}${metric}: got ${got}, want ${want}; "
    done
    if [[ -z "$bad" ]]; then
        record PASS "$name" "$# metrics as hand-counted"
    else
        record FAIL "$name" "$bad"
    fi
}

# --------------------------------------------------------------------------
# Cases: library, origin, called_cells
# --------------------------------------------------------------------------

OUT="$WORKDIR/stranded"
if run_metrics stranded "$OUT" --rrna-gtf "$FIX/added.gtf" --cell-barcodes "$FIX/barcodes.tsv" --strand Forward; then
    M="$OUT/metrics.csv"

    # 9 mapped reads in 12 mapped records: r4 and r5 add 3 secondary records.
    # r4's primary is on Mt_rRNA, so it is both an mtDNA and an rRNA read; r5's
    # secondary on chrM is not an mtDNA read, only an mtDNA alignment.
    expect library "$M" \
        "Mapped reads=9" \
        "Unmapped reads=1" \
        "Uniquely mapped reads=7" \
        "Multimapped reads (primary alignment)=2" \
        "Multimapped read alignments (all alignments)=5" \
        "rRNA reads in mapped reads (primary alignment)=3" \
        "Percentage of rRNA reads (of mapped reads, primary alignment)=0.3333" \
        "rRNA reads in uniquely mapped reads=2" \
        "Percentage of rRNA reads (of uniquely mapped reads)=0.2857" \
        "rRNA reads in multimapped reads (primary alignment)=1" \
        "Percentage of rRNA reads (of multimapped reads, primary alignment)=0.5000" \
        "rRNA alignments in multimapped read alignments (all alignments)=1" \
        "Percentage of rRNA alignments (of multimapped read alignments, all alignments)=0.2000" \
        "mtDNA reads in mapped reads (primary alignment)=4" \
        "Percentage of mtDNA reads (of mapped reads, primary alignment)=0.4444" \
        "mtDNA reads in uniquely mapped reads=3" \
        "Percentage of mtDNA reads (of uniquely mapped reads)=0.4286" \
        "mtDNA reads in multimapped reads (primary alignment)=1" \
        "Percentage of mtDNA reads (of multimapped reads, primary alignment)=0.5000" \
        "mtDNA alignments in multimapped read alignments (all alignments)=2" \
        "Percentage of mtDNA alignments (of multimapped read alignments, all alignments)=0.4000"

    # Of the 4 mtDNA reads: r1 and r4 sense, r2 antisense, r3 between genes
    expect origin "$M" \
        "mtDNA reads sense to mitochondrial genes (primary alignment)=2" \
        "Percentage of mtDNA reads sense to mitochondrial genes (of mtDNA reads)=0.5000" \
        "mtDNA reads antisense to mitochondrial genes (primary alignment)=1" \
        "Percentage of mtDNA reads antisense to mitochondrial genes (of mtDNA reads)=0.2500" \
        "mtDNA reads outside mitochondrial genes (primary alignment)=1" \
        "Percentage of mtDNA reads outside mitochondrial genes (of mtDNA reads)=0.2500"

    # AAAA holds r1 r2 r4 r5 r6 r10 (r8 is unmapped); TTTT is listed but absent
    expect called_cells "$M" \
        "Called cell barcodes=2" \
        "Called cell barcodes found in BAM=1" \
        "Mapped reads (called cells)=6" \
        "Uniquely mapped reads (called cells)=4" \
        "rRNA reads in mapped reads (primary alignment, called cells)=2" \
        "Percentage of rRNA reads (of mapped reads, primary alignment, called cells)=0.3333" \
        "rRNA reads in uniquely mapped reads (called cells)=1" \
        "Percentage of rRNA reads (of uniquely mapped reads, called cells)=0.2500" \
        "mtDNA reads in mapped reads (primary alignment, called cells)=3" \
        "Percentage of mtDNA reads (of mapped reads, primary alignment, called cells)=0.5000" \
        "mtDNA reads in uniquely mapped reads (called cells)=2" \
        "Percentage of mtDNA reads (of uniquely mapped reads, called cells)=0.5000"

    # Every row is a two-field CSV record although the names carry commas
    if awk -F'",' 'NR > 1 && NF != 2 { bad = 1 } END { exit bad }' "$M"; then
        record PASS "csv" "every row splits into name and value"
    else
        record FAIL "csv" "a row does not split into exactly name and value"
    fi

    # Library: sense 68+23 = 91, antisense 12+6 = 18 -> 18/109; exonic 12/80; intronic 6/29.
    # Called cells (AAAA; TTTT has no row): sense 80, antisense 15 -> 15/95.
    expect antisense "$OUT/antisense.csv" \
        "STARsolo strand=Forward" \
        "Reads sense to genes (exonic + intronic)=91" \
        "Reads antisense to genes (exonicAS + intronicAS)=18" \
        "Percentage of antisense reads (of reads assigned to genes)=0.1651" \
        "Percentage of exonic antisense reads (of exonic reads)=0.1500" \
        "Percentage of intronic antisense reads (of intronic reads)=0.2069" \
        "Called cell barcodes found in CellReads.stats=1" \
        "Reads sense to genes (exonic + intronic, called cells)=80" \
        "Reads antisense to genes (exonicAS + intronicAS, called cells)=15" \
        "Percentage of antisense reads (of reads assigned to genes, called cells)=0.1579" \
        "Percentage of exonic antisense reads (of exonic reads, called cells)=0.1429" \
        "Percentage of intronic antisense reads (of intronic reads, called cells)=0.2000"

    # Per barcode, primary alignments: AAAA r1 r2 r4 r5 r6 r10 (mtDNA r1 r2 r4, rRNA r4
    # on Mt_rRNA and r6), CCCC r3 r7 (mtDNA r3, rRNA r7 on the added contig), GGGG r9.
    # The called barcodes (AAAA) must add up to the called-cell rows above.
    if out="$("$PYTHON" - "$OUT/barcode_reads.tsv.gz" "$M" <<'PY'
import csv, gzip, sys
with gzip.open(sys.argv[1], "rt") as fh:
    rows = {r["CB"]: (int(r["MappedReads"]), int(r["MTReads"]), int(r["rRNAReads"]))
            for r in csv.DictReader(fh, delimiter="\t")}
want = {"AAAA": (6, 3, 2), "CCCC": (2, 1, 1), "GGGG": (1, 0, 0)}
bad = [f"{bc}: {rows.get(bc)} != {v}" for bc, v in want.items() if rows.get(bc) != v]
if set(rows) != set(want):
    bad.append(f"barcodes {sorted(rows)} != {sorted(want)}")
with open(sys.argv[2], newline="") as fh:
    m = {r[0]: r[1] for r in csv.reader(fh) if len(r) == 2}
cells = [rows.get("AAAA", (0, 0, 0))]
for i, key in enumerate(("Mapped reads (called cells)",
                         "mtDNA reads in mapped reads (primary alignment, called cells)",
                         "rRNA reads in mapped reads (primary alignment, called cells)")):
    if sum(c[i] for c in cells) != int(m[key]):
        bad.append(f"called-cell sum of column {i} != {key} ({m[key]})")
print("; ".join(bad) if bad else "ok")
sys.exit(1 if bad else 0)
PY
    )"; then
        record PASS "barcode_reads" "per-barcode reads as hand-counted; called barcodes sum to the called-cell rows"
    else
        record FAIL "barcode_reads" "$out"
    fi
fi

# --------------------------------------------------------------------------
# Cases: unstranded, no_cells
# --------------------------------------------------------------------------

OUT="$WORKDIR/unstranded"
if run_metrics unstranded "$OUT" --rrna-gtf "$FIX/added.gtf" --strand Unstranded; then
    M="$OUT/metrics.csv"
    expect unstranded "$M" \
        "mtDNA reads sense to mitochondrial genes (primary alignment)=N/A" \
        "mtDNA reads antisense to mitochondrial genes (primary alignment)=N/A" \
        "mtDNA reads outside mitochondrial genes (primary alignment)=1" \
        "Mapped reads=9"
    expect no_cells "$M" \
        "Called cell barcodes=N/A" \
        "Mapped reads (called cells)=N/A" \
        "Percentage of mtDNA reads (of mapped reads, primary alignment, called cells)=N/A"
    # STARsolo has no antisense direction to separate without a strand
    if [[ ! -e "$OUT/antisense.csv" ]]; then
        record PASS "antisense_unstranded" "no antisense file written"
    else
        record FAIL "antisense_unstranded" "antisense.csv written for an unstranded run"
    fi
fi

# --------------------------------------------------------------------------
# Case: percell
# --------------------------------------------------------------------------

# per-cell_images.py reads the stranded run's per-barcode counts, as PERCELL_METRICS
# reads CALC_READ_METRICS', so this case also checks that the two scripts fit together
COUNTS="$WORKDIR/stranded/barcode_reads.tsv.gz"
if ! "$PYTHON" -c 'import pandas, scipy, matplotlib' >/dev/null 2>&1; then
    record SKIP "percell" "a Python with pandas, scipy and matplotlib is needed"
elif [[ ! -s "$COUNTS" ]]; then
    record FAIL "percell" "no per-barcode counts from the stranded run: $COUNTS"
else
    SOLO="$FIX/sample_Solo.out"
    mkdir -p "$SOLO/GeneFull_Ex50pAS/raw" "$SOLO/Velocyto/raw" "$WORKDIR/percell"
    printf 'AAAA\nCCCC\nGGGG\n' > "$SOLO/GeneFull_Ex50pAS/raw/barcodes.tsv"
    # Intronic: AAAA 30 of 100, CCCC 1 of 10, GGGG none
    printf 'CB\tgenomeU\tgenomeM\texonic\tintronic\n' > "$SOLO/GeneFull_Ex50pAS/CellReads.stats"
    printf 'CBnotInPasslist\t10\t0\t3\t5\n' >> "$SOLO/GeneFull_Ex50pAS/CellReads.stats"
    printf 'AAAA\t100\t5\t60\t30\nCCCC\t10\t0\t5\t1\nGGGG\t0\t0\t0\t0\n' >> "$SOLO/GeneFull_Ex50pAS/CellReads.stats"
    # Velocyto, 2 genes x 3 barcodes: AAAA S6 U3 A1, CCCC U2, GGGG nothing
    printf 'AAAA\nCCCC\nGGGG\n' > "$SOLO/Velocyto/raw/barcodes.tsv"
    mtx() { printf '%%%%MatrixMarket matrix coordinate integer general\n%%\n2 3 %s\n' "$1"; shift; printf '%s\n' "$@"; }
    mtx 2 "1 1 4" "2 1 2" > "$SOLO/Velocyto/raw/spliced.mtx"
    mtx 2 "1 1 3" "2 2 2" > "$SOLO/Velocyto/raw/unspliced.mtx"
    mtx 1 "2 1 1"         > "$SOLO/Velocyto/raw/ambiguous.mtx"

    if run_logged "$BCA_TEST_LOGDIR/mt_rrna_metrics_percell.log" "$PYTHON" "$PERCELL_TOOL" \
            --solo-output "$SOLO" --read-counts "$COUNTS" --outdir "$WORKDIR/percell" \
            --cell-barcodes "$FIX/barcodes.tsv"; then
        # AAAA: 6 mapped, mtDNA r1 r2 r4, rRNA r4 (Mt_rRNA) r6. CCCC: r3 r7, r7 on the
        # added contig. GGGG: r9 only.
        if out="$("$PYTHON" - "$WORKDIR/percell/sample_metrics.csv" <<'PY'
import sys
import pandas as pd
df = pd.read_csv(sys.argv[1]).set_index("Cell")
want = {
    "AAAA": dict(MappedReads=6, MTPercent=50.0, rRNAPercent=100 * 2 / 6, IntronicPercent=30.0, UnsplicedPercent=30.0, IsCell=1),
    "CCCC": dict(MappedReads=2, MTPercent=50.0, rRNAPercent=50.0, IntronicPercent=10.0, UnsplicedPercent=100.0, IsCell=0),
    "GGGG": dict(MappedReads=1, MTPercent=0.0, rRNAPercent=0.0, IntronicPercent=0.0, UnsplicedPercent=0.0, IsCell=0),
}
bad = [f"{bc}.{col}: {df.loc[bc, col]} != {v}" for bc, cols in want.items()
       for col, v in cols.items() if abs(float(df.loc[bc, col]) - v) > 1e-6]
assert "TotalReads" not in df.columns, "legacy TotalReads column written"
print("; ".join(bad) if bad else "ok")
sys.exit(1 if bad else 0)
PY
        )"; then
            record PASS "percell" "per-cell reads, mtDNA, rRNA, intronic and unspliced as hand-counted"
        else
            record FAIL "percell" "$out"
        fi
    else
        record FAIL "percell" "per-cell_images.py failed, see mt_rrna_metrics_percell.log"
        tail_log "$BCA_TEST_LOGDIR/mt_rrna_metrics_percell.log"
    fi
fi

finish_check
