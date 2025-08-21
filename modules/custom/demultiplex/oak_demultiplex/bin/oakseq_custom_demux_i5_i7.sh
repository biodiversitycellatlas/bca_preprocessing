#!/bin/bash

# ------------------------------------------------------------------------------
# Usage / Argument parsing
# ------------------------------------------------------------------------------
barcode_raw=$1           # e.g. CAGGGTTGGC
index_type=$2            # either 'i5' or 'i7'
input_dir=$3
output_dir=$4
max_mm="${5:-1}"         # default: 1 mismatch

# Reverse‐complement function
revcomp() {
  echo "$1" | tr 'ACGTacgt' 'TGCAtgca' | rev
}

barcode_rc=$(revcomp "$barcode_raw")

# ------------------------------------------------------------------------------
# Paths & filenames
# ------------------------------------------------------------------------------
demux_dir="${output_dir}/fastq_demux_${index_type}_${barcode_raw}"
tmp_dir="${demux_dir}/tmp_demux_${index_type}_${barcode_raw}"

mkdir -p "${demux_dir}" "${tmp_dir}"

# Raw FASTQ files
R1="${input_dir}/Undetermined_S0_R1_001.fastq.gz"
R2="${input_dir}/Undetermined_S0_R2_001.fastq.gz"
I1="${input_dir}/Undetermined_S0_I1_001.fastq.gz"
I2="${input_dir}/Undetermined_S0_I2_001.fastq.gz"

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
zcat "${idx_fastq}" | \
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
  base=$(basename "${fq}" .fastq.gz)
  clean_base=${base#Undetermined_}
  seqtk subseq "${fq}" "${tmp_dir}/ids.txt" | gzip > "${demux_dir}/${barcode_raw}_${clean_base}.fastq.gz"
done

echo "Done."
echo "Outputs:"
ls -1 "${demux_dir}/${barcode_raw}"_*.fastq.gz
