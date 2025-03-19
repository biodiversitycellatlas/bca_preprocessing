#!/bin/bash

# ------------------------------------------------------------------
# Define Variables
# ------------------------------------------------------------------
work_dir=$1
mkdir -p ${work_dir}/Solo.out/Intronic/raw

exonic_mtx="${work_dir}/Solo.out/Gene/raw/matrix.mtx"
full_mtx="${work_dir}/Solo.out/GeneFull/raw/matrix.mtx"


# ------------------------------------------------------------------
# Capture the header lines and the numeric lines separately
# ------------------------------------------------------------------
# Get lines 1 and 2: MatrixMarket + comment line(s)
head -n 2 ${full_mtx} > gf_header1.txt

# Get the 3rd line: "num_features num_barcodes num_nonzero"
head -n 3 ${full_mtx} | tail -n 1 > gf_header2.txt

# Grab the coordinate lines (everything after the first 3 lines)
tail -n +4 ${full_mtx} > gf_body.txt
tail -n +4 ${exonic_mtx} > g_body.txt


# ------------------------------------------------------------------
# Calculate Intronic Counts
# Will paste them side by side and substract (GeneFull_count - Gene_count)
# Write the new 'count' to a new file intronic_body.txt
# ------------------------------------------------------------------
paste gf_body.txt g_body.txt \
| awk '{
    intronicCount = $3 - $6
    print $1, $2, intronicCount
}' > intronic_body.txt


# ------------------------------------------------------------------
# Create Intronic Matrix
# ------------------------------------------------------------------
# Recalculate how many nonzero lines are in intronic_body.txt
num_features=$(head -n 1 gf_header2.txt | awk '{print $1}')
num_barcodes=$(head -n 1 gf_header2.txt | awk '{print $2}')
num_nonzero=$(wc -l < intronic_body.txt)

# Write the matrix header
echo "%%MatrixMarket matrix coordinate integer general" > ${work_dir}/Solo.out/Intronic/raw/matrix.mtx
echo "%intronic count matrix derived from GeneFull - Gene" >> ${work_dir}/Solo.out/Intronic/raw/matrix.mtx
echo "${num_features} ${num_barcodes} ${num_nonzero}" >> ${work_dir}/Solo.out/Intronic/raw/matrix.mtx

# Append counts
cat intronic_body.txt >> ${work_dir}/Solo.out/Intronic/raw/matrix.mtx

# Copy features.tsv.gz and barcodes.tsv.gz so the structure lines up with the new matrix
cp ${work_dir}/Solo.out/Gene/raw/features.tsv ${work_dir}/Solo.out/Intronic/raw/features.tsv
cp ${work_dir}/Solo.out/Gene/raw/barcodes.tsv ${work_dir}/Solo.out/Intronic/raw/barcodes.tsv


# ------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------
rm gf_header1.txt gf_header2.txt gf_body.txt g_body.txt intronic_body.txt