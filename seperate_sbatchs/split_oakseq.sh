#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=05:00:00
#SBATCH --qos=short
#SBATCH --mem=32G
#SBATCH --job-name split_oakseq


# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
work_dir="/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec"
fastq_dir="${work_dir}/bcl2fastq2_nolanesplit_index"
output_dir="${work_dir}/fastq_splitted"
mkdir -p "${output_dir}"

# Defines i5 barcodes as an array
i5_barcodes=("AGGTAACACT" "AACGAACTGT" "TTTATCCTTG" "TTCGACAAGC")

# Set of FASTQs 
R1="${fastq_dir}/Undetermined_S0_R1_001.fastq.gz"
R2="${fastq_dir}/Undetermined_S0_R2_001.fastq.gz"
I2="${fastq_dir}/Undetermined_S0_I2_001.fastq.gz"

# Creates a temporary directory for read ID lists
tmp_dir="${work_dir}/tmp"
mkdir -p "${tmp_dir}"

# ------------------------------------------------------------------------------
# Extract read IDs from I2 allowing for 1 mismatch, saving non-matches separately
# ------------------------------------------------------------------------------
echo "Extracting read IDs from I2 with up to 1 mismatch allowed..."
zcat "${I2}" | awk -v bc1="${i5_barcodes[0]}" -v bc2="${i5_barcodes[1]}" \
                         -v bc3="${i5_barcodes[2]}" -v bc4="${i5_barcodes[3]}" \
                         -v tmp_dir="${tmp_dir}" '
function hamming(a, b) {
    if (length(a) != length(b)) return -1;
    d = 0;
    for (i = 1; i <= length(a); i++) {
        if (substr(a, i, 1) != substr(b, i, 1))
            d++;
    }
    return d;
}
NR % 4 == 1 {
    # Extract read ID (remove the "@" and any trailing extra info)
    header = $0; sub(/^@/, "", header); split(header, a, " ");
    readid = a[1]
}
NR % 4 == 2 {
    matched = "";
    if (hamming($0, bc1) <= 1) matched = bc1;
    else if (hamming($0, bc2) <= 1) matched = bc2;
    else if (hamming($0, bc3) <= 1) matched = bc3;
    else if (hamming($0, bc4) <= 1) matched = bc4;
    if (matched != "") {
      print readid > tmp_dir"/"matched".ids.txt";
    } else {
      print readid > tmp_dir"/unfiltered.ids.txt";
    }
}
'

# ------------------------------------------------------------------------------
# Extract corresponding reads using seqtk subseq
# ------------------------------------------------------------------------------
echo "Extracting FASTQ subsets with seqtk..."
# Loop over each barcode and the unfiltered category
for bc in "${i5_barcodes[@]}" "unfiltered"; do
    id_file="${tmp_dir}/${bc}.ids.txt"
    if [[ -s "${id_file}" ]]; then
        seqtk subseq "${R1}" "${id_file}" | gzip > "${output_dir}/${bc}_R1.fastq.gz"
        seqtk subseq "${R2}" "${id_file}" | gzip > "${output_dir}/${bc}_R2.fastq.gz"
        seqtk subseq "${I2}" "${id_file}" | gzip > "${output_dir}/${bc}_I2.fastq.gz"
        echo "Extracted reads for ${bc}"
    else
        echo "No reads found for ${bc}"
    fi
done

# Clean up temporary files
rm -rf "${tmp_dir}"

echo "Done. Split FASTQs are in: ${output_dir}"