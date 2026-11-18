#!/usr/bin/env bash
# Split-barcode variant of the seqspec demo, for an end-to-end run of a CB_UMI_Complex layout.
#
# Rewrites read 1 of the 10x example (16-nt barcode + 12-nt UMI) as
#     barcode[1-8] + linker TTGCAC + barcode[9-16] + UMI        (34 nt)
# and writes a seqspec file with two barcode regions whose local onlists are the distinct
# first and second halves of the 10x v3 whitelist. Joining a read's two halves gives back its
# 10x barcode, so the run can be compared with the single-barcode demo.
#
# Usage:
#   make_seqspec_split_fixture.sh --repo REPO_ROOT --whitelist 3M-february-2018.txt[.gz] --outdir DIR
#
# Writes into DIR: split_R1.fastq.gz, split_R2.fastq.gz, bc_half1.txt, bc_half2.txt,
# split.seqspec.yaml, samplesheet.csv and params.config (the run's -c config).
set -euo pipefail

repo="" whitelist="" outdir=""
while (( $# )); do
    case "$1" in
        --repo)      repo="$2"; shift 2 ;;
        --whitelist) whitelist="$2"; shift 2 ;;
        --outdir)    outdir="$2"; shift 2 ;;
        -h|--help)   sed -n '2,15p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$repo" && -n "$whitelist" && -n "$outdir" ]] || { echo "need --repo, --whitelist and --outdir" >&2; exit 2; }

data="$repo/examples/10x/test_data"
linker="TTGCAC"
mkdir -p "$outdir"
outdir=$(cd "$outdir" && pwd)

# Read 1: insert the linker after the 8th barcode base (quality 'I' for its bases)
gzip -dc "$data/spis_10x_R1.fastq.gz" \
    | awk -v L="$linker" 'NR % 4 == 2 { $0 = substr($0, 1, 8) L substr($0, 9) }
                          NR % 4 == 0 { $0 = substr($0, 1, 8) "IIIIII" substr($0, 9) }
                          NR % 4 == 1 { sub(/length=[0-9]+/, "length=34") }
                          { print }' \
    | gzip -c > "$outdir/split_R1.fastq.gz"
cp "$data/spis_10x_R2.fastq.gz" "$outdir/split_R2.fastq.gz"

# Onlists: the distinct halves of the whitelist
gzip -dcf "$whitelist" | tr -d '\r' | awk 'NF { print substr($1, 1, 8) }' | sort -u > "$outdir/bc_half1.txt"
gzip -dcf "$whitelist" | tr -d '\r' | awk 'NF { print substr($1, 9, 8) }' | sort -u > "$outdir/bc_half2.txt"

md5() { md5sum "$1" | cut -d' ' -f1; }
size() { wc -c < "$1" | tr -d ' '; }

cat > "$outdir/split.seqspec.yaml" <<EOF
!Assay
seqspec_version: 0.4.0
assay_id: spis-split-barcode-test
name: Split-barcode test layout (10x example reads, barcode split by a linker)
doi: ''
date: 18 November 2026
description: >-
  Test fixture written by tests/lib/make_seqspec_split_fixture.sh. Read 1 of the 10x
  example with a 6-nt linker inserted in the middle of the cell barcode.
modalities:
- rna
lib_struct: ''
sequence_protocol: Illumina NextSeq 500 (EFO:0009173)
sequence_kit: NextSeq 500/550 High-Output v2.5 Kit
library_protocol: single-cell RNA sequencing assay (OBI:0002631)
library_kit: test fixture
sequence_spec:
- !Read
  read_id: R1
  name: Read 1
  modality: rna
  primer_id: truseq_read1
  min_len: 34
  max_len: 34
  strand: pos
  files:
  - !File
    file_id: split_R1.fastq.gz
    filename: split_R1.fastq.gz
    filetype: fastq
    filesize: $(size "$outdir/split_R1.fastq.gz")
    url: split_R1.fastq.gz
    urltype: local
    md5: $(md5 "$outdir/split_R1.fastq.gz")
- !Read
  read_id: R2
  name: Read 2
  modality: rna
  primer_id: truseq_read2
  min_len: 56
  max_len: 56
  strand: neg
  files:
  - !File
    file_id: split_R2.fastq.gz
    filename: split_R2.fastq.gz
    filetype: fastq
    filesize: $(size "$outdir/split_R2.fastq.gz")
    url: split_R2.fastq.gz
    urltype: local
    md5: $(md5 "$outdir/split_R2.fastq.gz")
library_spec:
- !Region
  region_id: rna
  region_type: rna
  name: rna
  sequence_type: joined
  sequence: AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNN${linker}NNNNNNNNXXXXXXXXXXXXXAGATCGGAAGAGCACACGTCTGAACTCCAGTCACXXXXXXXXATCTCGTATGCCGTCTTCTGCTTG
  min_len: 163
  max_len: 260
  onlist: null
  regions:
  - !Region
    region_id: illumina_p5
    region_type: illumina_p5
    name: Illumina P5
    sequence_type: fixed
    sequence: AATGATACGGCGACCACCGAGATCTACAC
    min_len: 29
    max_len: 29
    onlist: null
    regions: []
  - !Region
    region_id: truseq_read1
    region_type: truseq_read1
    name: TruSeq Read 1
    sequence_type: fixed
    sequence: ACACTCTTTCCCTACACGACGCTCTTCCGATCT
    min_len: 33
    max_len: 33
    onlist: null
    regions: []
  - !Region
    region_id: barcode_1
    region_type: barcode
    name: Cell barcode, first half
    sequence_type: onlist
    sequence: NNNNNNNN
    min_len: 8
    max_len: 8
    onlist: !Onlist
      file_id: bc_half1.txt
      filename: bc_half1.txt
      filetype: txt
      filesize: $(size "$outdir/bc_half1.txt")
      url: bc_half1.txt
      urltype: local
      md5: $(md5 "$outdir/bc_half1.txt")
    regions: []
  - !Region
    region_id: linker
    region_type: linker
    name: Linker
    sequence_type: fixed
    sequence: ${linker}
    min_len: 6
    max_len: 6
    onlist: null
    regions: []
  - !Region
    region_id: barcode_2
    region_type: barcode
    name: Cell barcode, second half
    sequence_type: onlist
    sequence: NNNNNNNN
    min_len: 8
    max_len: 8
    onlist: !Onlist
      file_id: bc_half2.txt
      filename: bc_half2.txt
      filetype: txt
      filesize: $(size "$outdir/bc_half2.txt")
      url: bc_half2.txt
      urltype: local
      md5: $(md5 "$outdir/bc_half2.txt")
    regions: []
  - !Region
    region_id: umi
    region_type: umi
    name: UMI
    sequence_type: random
    sequence: XXXXXXXXXXXX
    min_len: 12
    max_len: 12
    onlist: null
    regions: []
  - !Region
    region_id: cdna
    region_type: cdna
    name: cDNA
    sequence_type: random
    sequence: X
    min_len: 1
    max_len: 98
    onlist: null
    regions: []
  - !Region
    region_id: truseq_read2
    region_type: truseq_read2
    name: TruSeq Read 2
    sequence_type: fixed
    sequence: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    min_len: 34
    max_len: 34
    onlist: null
    regions: []
  - !Region
    region_id: index7
    region_type: index7
    name: i7 index
    sequence_type: random
    sequence: XXXXXXXX
    min_len: 8
    max_len: 8
    onlist: null
    regions: []
  - !Region
    region_id: illumina_p7
    region_type: illumina_p7
    name: Illumina P7
    sequence_type: fixed
    sequence: ATCTCGTATGCCGTCTTCTGCTTG
    min_len: 24
    max_len: 24
    onlist: null
    regions: []
EOF

printf 'sample,fastq_cDNA,fastq_BC_UMI,expected_cells,p5,p7,rt\nspis_split,split_R2.fastq.gz,split_R1.fastq.gz,30,,,\n' > "$outdir/samplesheet.csv"

cat > "$outdir/params.config" <<EOF
params {
    input = "$outdir/samplesheet.csv"
    outdir = "$outdir/output/"
    ref_gtf = "\${projectDir}/examples/10x/test_data/spis_long.annot.gtf"
    ref_fasta = "\${projectDir}/examples/10x/test_data/spis_gDNA.fa"
    protocol = "seqspec"
    seqspec_file = "$outdir/split.seqspec.yaml"
    mapping_software = "both"
    perform_geneext = false
    perform_10x_saturate = true
    perform_doublet_detection = false
    perform_kraken = false
}
EOF

echo "Split-barcode fixture in $outdir: $(wc -l < "$outdir/bc_half1.txt") + $(wc -l < "$outdir/bc_half2.txt") onlist entries"
