#!/bin/bash

# ------------------------------------------------------------------------------
# Usage / Argument parsing
# ------------------------------------------------------------------------------
barcode_raw=$1           # e.g. CAGGGTTGGC
index_type=$2            # either 'i5' or 'i7'
input_dir=$3
outdir=$4
max_mm="${5:-1}"         # default: 1 mismatch

# Reverse‐complement function
revcomp() {
  echo "$1" | tr 'ACGTacgt' 'TGCAtgca' | rev
}

barcode_rc=$(revcomp "$barcode_raw")

# ------------------------------------------------------------------------------
# Paths & filenames
# ------------------------------------------------------------------------------
demux_dir="${outdir}/fastq_demux_${index_type}_${barcode_raw}"
tmp_dir="${demux_dir}/tmp_demux_${index_type}_${barcode_raw}"

mkdir -p "${demux_dir}" "${tmp_dir}"

# Raw FASTQ files: accept plain or gzipped, whichever is present
find_read() {
  local read_id=$1
  local candidate
  for ext in fastq.gz fq.gz fastq fq; do
    candidate="${input_dir}/Undetermined_S0_${read_id}_001.${ext}"
    if [[ -f "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done
  echo "Error: no FASTQ found for ${read_id} in ${input_dir}" >&2
  return 1
}

R1=$(find_read R1) || exit 1
R2=$(find_read R2) || exit 1
I1=$(find_read I1) || exit 1
I2=$(find_read I2) || exit 1

# Stream a FASTQ as plain text, whether or not it is gzipped
cat_fastq() {
  if [[ "$(od -An -N2 -tx1 "$1" | tr -d ' \n')" == "1f8b" ]]; then
    gzip -dc "$1"
  else
    cat "$1"
  fi
}

# Pick the index FASTQ based on i5 or i7
if [[ "${index_type}" == "i5" ]]; then
  idx_fastq="${I2}"
elif [[ "${index_type}" == "i7" ]]; then
  idx_fastq="${I1}"
else
  echo "Error: first argument must be 'i5' or 'i7'."
  exit 1
fi

# ------------------------------------------------------------------------------
# 1) Extract matching read IDs
# ------------------------------------------------------------------------------
echo "Extracting read IDs from ${index_type} (${idx_fastq}) for barcode ${barcode_raw} (RC=${barcode_rc}), ≤${max_mm} mismatches..."
cat_fastq "${idx_fastq}" | \
awk -v bc="${barcode_rc}" -v mm="${max_mm}" -v out="${tmp_dir}/ids.txt" '
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
    if (hamming($0,bc) >= 0 && hamming($0,bc) <= mm) {
      print id >> out
    }
  }
'

# ------------------------------------------------------------------------------
# 2) Subset all four FASTQs with seqtk
# ------------------------------------------------------------------------------
echo "Demultiplexing all four reads into ${demux_dir}/${barcode_raw}_*.fastq.gz …"
for fq in "${R1}" "${R2}" "${I1}" "${I2}"; do
  # seqtk reads plain and gzipped FASTQs alike; strip whichever extension is present
  base=$(basename "${fq}" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
  clean_base=${base#Undetermined_}
  seqtk subseq "${fq}" "${tmp_dir}/ids.txt" | gzip > "${demux_dir}/${barcode_raw}_${clean_base}.fastq.gz"
done

echo "Done."
echo "Outputs:"
ls -1 "${demux_dir}/${barcode_raw}"_*.fastq.gz
