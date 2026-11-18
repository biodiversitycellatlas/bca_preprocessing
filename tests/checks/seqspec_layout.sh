#!/usr/bin/env bash
# description: Check the STARsolo and alevin-fry barcode layouts bin/seqspec_layout.py derives from seqspec files
#
# Protocol 'seqspec' reads the barcode and UMI positions of both mappers from the
# seqspec file (SEQSPEC_LAYOUT runs bin/seqspec_layout.py once per run). The script
# replaces seqspec's own 'index -t starsolo', which
#   - prints --soloType itself, so the '--soloType CB_UMI_Simple $(seqspec index ...)'
#     the pipeline used gave STAR the parameter twice (fatal: duplicate parameter);
#   - writes --soloCBlen as the barcode's end, wrong when the barcode does not start the read;
#   - keeps only the first barcode region, so split barcodes were impossible.
#
# Each case writes a small seqspec file, runs the script and compares the layout,
# the alevin geometry and the onlist table with values worked out by hand.
#
# Needs a Python with the seqspec package (the STARSOLO_ALIGN conda env has it);
# no sequencing data, STAR or Nextflow.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

TOOL="$PROJECT_ROOT/bin/seqspec_layout.py"
DEMO_SPEC="$PROJECT_ROOT/examples/seqspec/conf/spis_10xv3.seqspec.yaml"
PYTHON="${BCA_PYTHON:-}"
KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh seqspec_layout [-- OPTIONS]
       tests/checks/seqspec_layout.sh [OPTIONS]

Validate the barcode layouts derived from seqspec files.

Options:
  --python PATH   Python interpreter with the seqspec package (default: \$BCA_PYTHON,
                  else python3). Its folder is put first on PATH, for 'seqspec check'.
  --keep          Keep the generated specs and outputs for inspection.
  -h, --help      Show this message.

Cases:
  demo_10xv3      the demo spec (examples/seqspec): CB_UMI_Simple 1/16/17/12,
                  1[1-16] / 1[17-28] / 2[1-end], one remote onlist
  spacer          a 4-nt spacer before the barcode and before the cDNA:
                  --soloCBstart 5 --soloCBlen 16, 1[5-20], 2[5-end]
  split           three barcodes between linkers: CB_UMI_Complex positions,
                  soloCBmatchWLtype 1MM, 1[1-9,14-22,27-35], onlists in read order
  env_roundtrip   the env file sources in bash and gives back the values
  variable_len    a variable-length region before the barcode is refused
  umi_other_read  barcode and UMI in different reads are refused
  no_umi          a layout without UMI is refused
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

log_header "seqspec_layout"

if [[ -z "$PYTHON" ]]; then
    PYTHON=$(command -v python3 || command -v python || true)
fi
if [[ -z "$PYTHON" ]] || ! "$PYTHON" -c 'import seqspec' 2>/dev/null; then
    record SKIP "all cases" "no Python with the seqspec package (pass --python, e.g. the STARSOLO_ALIGN conda env's)"
    finish_check
fi
export PATH="$(dirname "$PYTHON"):$PATH"

WORKDIR="$(mktemp -d -t bca_seqspec_layout.XXXXXX)"
if (( KEEP )); then
    log_dim "keeping $WORKDIR"
else
    trap 'rm -rf "$WORKDIR"' EXIT
fi

# --------------------------------------------------------------------------
# Spec writer
# --------------------------------------------------------------------------

# write_spec FILE R1_LEN R2_LEN REGION...
#   REGION is id:type:sequence_type:sequence:min:max[:onlist_url:urltype], in library
#   order (5' to 3'). Read 1 starts after the region with id truseq_read1 (strand pos),
#   read 2 reads back from the region with id truseq_read2 (strand neg).
write_spec() {
    local out="$1" r1_len="$2" r2_len="$3"; shift 3
    local joined="" min_sum=0 max_sum=0 regions="" def
    for def in "$@"; do
        # The onlist URL may hold ':' itself (https://...), so the urltype is split off the end
        IFS=: read -r id type stype seq min max tail <<<"$def"
        url="" urltype=""
        if [[ -n "${tail:-}" ]]; then url="${tail%:*}"; urltype="${tail##*:}"; fi
        joined+="$seq"; min_sum=$(( min_sum + min )); max_sum=$(( max_sum + max ))
        regions+="  - !Region
    region_id: $id
    region_type: $type
    name: $id
    sequence_type: $stype
    sequence: $seq
    min_len: $min
    max_len: $max
"
        if [[ -n "$url" ]]; then
            regions+="    onlist: !Onlist
      file_id: $id.onlist
      filename: $(basename "$url")
      filetype: txt
      filesize: 0
      url: $url
      urltype: $urltype
      md5: d41d8cd98f00b204e9800998ecf8427e
"
        else
            regions+="    onlist: null
"
        fi
        regions+="    regions: []
"
    done
    cat > "$out" <<EOF
!Assay
seqspec_version: 0.4.0
assay_id: test
name: test
doi: ''
date: 18 November 2026
description: seqspec_layout check fixture
modalities:
- rna
lib_struct: ''
sequence_protocol: Illumina NextSeq 500 (EFO:0009173)
sequence_kit: NextSeq 500/550 High-Output v2.5 Kit
library_protocol: single-cell RNA sequencing assay (OBI:0002631)
library_kit: test
sequence_spec:
- !Read
  read_id: R1
  name: Read 1
  modality: rna
  primer_id: truseq_read1
  min_len: $r1_len
  max_len: $r1_len
  strand: pos
  files: []
- !Read
  read_id: R2
  name: Read 2
  modality: rna
  primer_id: truseq_read2
  min_len: $r2_len
  max_len: $r2_len
  strand: neg
  files: []
library_spec:
- !Region
  region_id: rna
  region_type: rna
  name: rna
  sequence_type: joined
  sequence: $joined
  min_len: $min_sum
  max_len: $max_sum
  onlist: null
  regions:
$regions
EOF
}

P5="illumina_p5:illumina_p5:fixed:AATGATACGGCGACCACCGAGATCTACAC:29:29"
R1P="truseq_read1:truseq_read1:fixed:ACACTCTTTCCCTACACGACGCTCTTCCGATCT:33:33"
CDNA="cdna:cdna:random:X:1:98"
R2P="truseq_read2:truseq_read2:fixed:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC:34:34"
P7="illumina_p7:illumina_p7:fixed:ATCTCGTATGCCGTCTTCTGCTTG:24:24"

# run_layout NAME SPEC -- runs the tool in WORKDIR/NAME; the exit status is the tool's
run_layout() {
    local name="$1" spec="$2"
    mkdir -p "$WORKDIR/$name"
    (cd "$WORKDIR/$name" && "$PYTHON" "$TOOL" --spec "$spec" --no-check) >"$WORKDIR/$name/log" 2>&1
}

# value NAME KEY -- one setting of a case's env file, as bash reads it
value() {
    (set +u; source "$WORKDIR/$1/seqspec_layout.env"; eval "printf '%s' \"\${$2}\"")
}

# expect NAME KEY WANTED -- empty output when equal, else a mismatch note
expect() {
    local got
    got=$(value "$1" "$2")
    [[ "$got" == "$3" ]] || printf '%s: got "%s", expected "%s"; ' "$2" "$got" "$3"
}

# expect_fail NAME SPEC PATTERN -- the tool must fail with a message matching PATTERN
expect_fail() {
    local name="$1" spec="$2" pattern="$3"
    if run_layout "$name" "$spec"; then
        record FAIL "$name" "accepted the spec"
    elif grep -q -- "$pattern" "$WORKDIR/$name/log"; then
        record PASS "$name" "$(grep -m1 '^ERROR' "$WORKDIR/$name/log" | cut -c1-110)"
    else
        record FAIL "$name" "failed without '$pattern': $(tail -n1 "$WORKDIR/$name/log")"
    fi
}

# --------------------------------------------------------------------------
# Cases
# --------------------------------------------------------------------------

# demo_10xv3: the spec the demo runs with
if [[ ! -f "$DEMO_SPEC" ]]; then
    record FAIL "demo_10xv3" "$DEMO_SPEC is missing"
elif ! run_layout demo_10xv3 "$DEMO_SPEC"; then
    record FAIL "demo_10xv3" "$(tail -n1 "$WORKDIR/demo_10xv3/log")"
else
    msg=$(expect demo_10xv3 SEQSPEC_STAR_TYPE "--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0")
    msg+=$(expect demo_10xv3 SEQSPEC_STAR_CBMATCH "")
    msg+=$(expect demo_10xv3 SEQSPEC_ALEVIN_BC "1[1-16]")
    msg+=$(expect demo_10xv3 SEQSPEC_ALEVIN_UMI "1[17-28]")
    msg+=$(expect demo_10xv3 SEQSPEC_ALEVIN_READ "2[1-end]")
    msg+=$(expect demo_10xv3 SEQSPEC_BARCODE_READ "R1")
    msg+=$(expect demo_10xv3 SEQSPEC_CDNA_READ "R2")
    onlists=$(cat "$WORKDIR/demo_10xv3/onlists.tsv")
    [[ "$onlists" == $'1\tbarcode\thttps\thttps://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz' ]] \
        || msg+="onlists.tsv: '$onlists'; "
    if [[ -z "$msg" ]]; then record PASS "demo_10xv3"; else record FAIL "demo_10xv3" "$msg"; fi
fi

# spacer: barcode after a fixed 4-nt spacer, cDNA read starting with a 4-nt linker
write_spec "$WORKDIR/spacer.yaml" 32 60 "$P5" "$R1P" \
    "spacer:linker:fixed:ACGT:4:4" "barcode:barcode:onlist:NNNNNNNNNNNNNNNN:16:16" \
    "umi:umi:random:XXXXXXXXXXXX:12:12" "$CDNA" "tso:linker:fixed:TTTT:4:4" "$R2P" "$P7"
if ! run_layout spacer "$WORKDIR/spacer.yaml"; then
    record FAIL "spacer" "$(tail -n1 "$WORKDIR/spacer/log")"
else
    msg=$(expect spacer SEQSPEC_STAR_TYPE "--soloType CB_UMI_Simple --soloCBstart 5 --soloCBlen 16 --soloUMIstart 21 --soloUMIlen 12 --soloBarcodeReadLength 0")
    msg+=$(expect spacer SEQSPEC_ALEVIN_BC "1[5-20]")
    msg+=$(expect spacer SEQSPEC_ALEVIN_UMI "1[21-32]")
    msg+=$(expect spacer SEQSPEC_ALEVIN_READ "2[5-end]")
    # The barcode region has no onlist: one row with empty urltype and url
    [[ "$(cat "$WORKDIR/spacer/onlists.tsv")" == $'1\tbarcode\t\t' ]] || msg+="onlists.tsv: '$(cat "$WORKDIR/spacer/onlists.tsv")'; "
    seqspec_own=$(seqspec index -m rna -t starsolo -s read -i R1 "$WORKDIR/spacer.yaml" 2>/dev/null || true)
    if [[ -z "$msg" ]]; then record PASS "spacer" "seqspec's own: ${seqspec_own:-n/a}"; else record FAIL "spacer" "$msg"; fi
fi

# split: three 9-nt barcodes separated by 4-nt linkers, then an 8-nt UMI (43 nt)
write_spec "$WORKDIR/split.yaml" 43 56 "$P5" "$R1P" \
    "bc1:barcode:onlist:NNNNNNNNN:9:9:lists/bc1.txt:local" "link1:linker:fixed:ACGT:4:4" \
    "bc2:barcode:onlist:NNNNNNNNN:9:9:lists/bc2.txt:local" "link2:linker:fixed:TGCA:4:4" \
    "bc3:barcode:onlist:NNNNNNNNN:9:9:https://example.org/bc3.txt.gz:https" \
    "umi:umi:random:XXXXXXXX:8:8" "$CDNA" "$R2P" "$P7"
if ! run_layout split "$WORKDIR/split.yaml"; then
    record FAIL "split" "$(tail -n1 "$WORKDIR/split/log")"
else
    msg=$(expect split SEQSPEC_STAR_TYPE "--soloType CB_UMI_Complex --soloCBposition 0_0_0_8 0_13_0_21 0_26_0_34 --soloUMIposition 0_35_0_42")
    msg+=$(expect split SEQSPEC_STAR_CBMATCH "1MM")
    msg+=$(expect split SEQSPEC_ALEVIN_BC "1[1-9,14-22,27-35]")
    msg+=$(expect split SEQSPEC_ALEVIN_UMI "1[36-43]")
    msg+=$(expect split SEQSPEC_N_BARCODES "3")
    want=$'1\tbc1\tlocal\tlists/bc1.txt\n2\tbc2\tlocal\tlists/bc2.txt\n3\tbc3\thttps\thttps://example.org/bc3.txt.gz'
    [[ "$(cat "$WORKDIR/split/onlists.tsv")" == "$want" ]] || msg+="onlists.tsv: '$(tr '\n\t' '|,' < "$WORKDIR/split/onlists.tsv")'; "
    [[ ! -e "$WORKDIR/split/onlist_3.txt" ]] || msg+="downloaded an onlist without --fetch-onlists; "
    if [[ -z "$msg" ]]; then record PASS "split"; else record FAIL "split" "$msg"; fi
fi

# env_roundtrip: the values survive shell quoting (brackets, commas, spaces)
if [[ -f "$WORKDIR/split/seqspec_layout.env" ]]; then
    got=$(bash -c 'set -eu; source "$1"; printf "%s|%s" "$SEQSPEC_STAR_TYPE" "$SEQSPEC_ALEVIN_BC"' _ "$WORKDIR/split/seqspec_layout.env" 2>&1 || true)
    if [[ "$got" == "--soloType CB_UMI_Complex --soloCBposition 0_0_0_8 0_13_0_21 0_26_0_34 --soloUMIposition 0_35_0_42|1[1-9,14-22,27-35]" ]]; then
        record PASS "env_roundtrip"
    else
        record FAIL "env_roundtrip" "$got"
    fi
else
    record SKIP "env_roundtrip" "the split case wrote no env file"
fi

# variable_len: a 0-4 nt spacer before the barcode leaves no fixed positions
write_spec "$WORKDIR/variable.yaml" 32 56 "$P5" "$R1P" \
    "stagger:linker:random:XXXX:0:4" "barcode:barcode:onlist:NNNNNNNNNNNNNNNN:16:16" \
    "umi:umi:random:XXXXXXXXXXXX:12:12" "$CDNA" "$R2P" "$P7"
expect_fail variable_len "$WORKDIR/variable.yaml" "variable length"

# umi_other_read: barcode in read 1, UMI read first by read 2
write_spec "$WORKDIR/umi_r2.yaml" 16 56 "$P5" "$R1P" \
    "barcode:barcode:onlist:NNNNNNNNNNNNNNNN:16:16" "$CDNA" "umi:umi:random:XXXXXXXXXXXX:12:12" "$R2P" "$P7"
expect_fail umi_other_read "$WORKDIR/umi_r2.yaml" "both must be in the same read"

# no_umi
write_spec "$WORKDIR/no_umi.yaml" 16 56 "$P5" "$R1P" \
    "barcode:barcode:onlist:NNNNNNNNNNNNNNNN:16:16" "$CDNA" "$R2P" "$P7"
expect_fail no_umi "$WORKDIR/no_umi.yaml" "covers a 'umi' region"

finish_check
