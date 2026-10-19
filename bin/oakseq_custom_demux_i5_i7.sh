#!/bin/bash

# ------------------------------------------------------------------------------
# Demultiplex OAK reads on one index (i5 or i7)
#
# Usage:
#   oakseq_custom_demux_i5_i7.sh --barcode CAGGGTTGGC --index-type i5|i7 \
#       --r1 R1.fastq.gz --r2 R2.fastq.gz [--i1 I1.fastq.gz] --i2 I2.fastq.gz \
#       --out PREFIX [--max-mm 1]
#
# Reads whose index (I2 for i5, I1 for i7) lies within --max-mm mismatches of the
# reverse-complemented barcode are written to PREFIX_{R1,R2,I1,I2}_001.fastq.gz.
# --i1 is only required when demultiplexing on i7.
# ------------------------------------------------------------------------------
set -euo pipefail

usage() {
  sed -n '4,13p' "$0" | sed 's/^# \{0,1\}//'
  exit 1
}

# ------------------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------------------
barcode_raw=""
index_type=""
R1=""
R2=""
I1=""
I2=""
out_prefix=""
max_mm=1                 # default: 1 mismatch

while [[ $# -gt 0 ]]; do
  case "$1" in
    --barcode)    barcode_raw=$2; shift 2 ;;
    --index-type) index_type=$2;  shift 2 ;;
    --r1)         R1=$2;          shift 2 ;;
    --r2)         R2=$2;          shift 2 ;;
    --i1)         I1=$2;          shift 2 ;;
    --i2)         I2=$2;          shift 2 ;;
    --out)        out_prefix=$2;  shift 2 ;;
    --max-mm)     max_mm=$2;      shift 2 ;;
    -h|--help)    usage ;;
    *) echo "Error: unknown argument '$1'" >&2; usage ;;
  esac
done

if [[ -z "${barcode_raw}" || -z "${index_type}" || -z "${R1}" || -z "${R2}" || -z "${out_prefix}" ]]; then
  echo "Error: --barcode, --index-type, --r1, --r2 and --out are required." >&2
  usage
fi

# Pick the index FASTQ based on i5 or i7
if [[ "${index_type}" == "i5" ]]; then
  idx_fastq="${I2}"
elif [[ "${index_type}" == "i7" ]]; then
  idx_fastq="${I1}"
else
  echo "Error: --index-type must be 'i5' or 'i7'." >&2
  exit 1
fi

if [[ -z "${idx_fastq}" ]]; then
  echo "Error: demultiplexing on ${index_type} needs --$([[ ${index_type} == i5 ]] && echo i2 || echo i1)." >&2
  exit 1
fi

for fq in "${R1}" "${R2}" "${I1}" "${I2}"; do
  if [[ -n "${fq}" && ! -f "${fq}" ]]; then
    echo "Error: FASTQ '${fq}' not found." >&2
    exit 1
  fi
done

# Reverse‐complement function
revcomp() {
  echo "$1" | tr 'ACGTacgt' 'TGCAtgca' | rev
}

barcode_rc=$(revcomp "$barcode_raw")

mkdir -p "$(dirname "${out_prefix}")"
ids_file="${out_prefix}_${index_type}_ids.txt"
: > "${ids_file}"

# Stream a FASTQ as plain text, whether or not it is gzipped
cat_fastq() {
  if [[ "$(od -An -N2 -tx1 "$1" | tr -d ' \n')" == "1f8b" ]]; then
    gzip -dc "$1"
  else
    cat "$1"
  fi
}

# ------------------------------------------------------------------------------
# 1) Extract matching read IDs
# ------------------------------------------------------------------------------
echo "Extracting read IDs from ${index_type} (${idx_fastq}) for barcode ${barcode_raw} (RC=${barcode_rc}), ≤${max_mm} mismatches..."
cat_fastq "${idx_fastq}" | \
awk -v bc="${barcode_rc}" -v mm="${max_mm}" -v out="${ids_file}" '
  function hamming(a,b) {
    if (length(a)!=length(b)) return -1;
    d=0;
    for(i=1;i<=length(a);i++) if(substr(a,i,1)!=substr(b,i,1)) d++;
    return d;
  }
  NR%4==1 {
    hdr=$0; sub(/^@/,"",hdr); split(hdr,A," "); id=A[1];
  }
  NR%4==2 {
    h=hamming($0,bc);
    if (h >= 0 && h <= mm) {
      print id >> out
    }
  }
'

n_matched=$(wc -l < "${ids_file}" | tr -d ' ')
echo "Reads matching ${index_type} ${barcode_raw}: ${n_matched}"

# ------------------------------------------------------------------------------
# 2) Subset all FASTQs with seqtk
# ------------------------------------------------------------------------------
echo "Demultiplexing reads into ${out_prefix}_*_001.fastq.gz …"
subset() {
  local fq=$1
  local read_id=$2
  # seqtk reads plain and gzipped FASTQs alike
  seqtk subseq "${fq}" "${ids_file}" | gzip > "${out_prefix}_${read_id}_001.fastq.gz"
}

subset "${R1}" R1
subset "${R2}" R2
if [[ -n "${I1}" ]]; then subset "${I1}" I1; fi
if [[ -n "${I2}" ]]; then subset "${I2}" I2; fi

rm -f "${ids_file}"

echo "Done."
echo "Outputs:"
ls -1 "${out_prefix}"_*_001.fastq.gz
