#!/usr/bin/env python

import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import seaborn as sns


# ------------------------------------------------------------------
# LOAD FILES
# ------------------------------------------------------------------
# dir_list = ["/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002/mapping_STARsolo/mapping_STARsolo_N/", 
#                   "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA003_BCA004_FINAL/mapping_STARsolo/mapping_STARsolo_N/",
#                   "/users/asebe/bvanwaardenburg/git/data/241204_ParseBio_Nvec_Tcas_nuclei/Nvec_FINAL/mapping_STARsolo/mapping_STARsolo_N/"]

# dir_list = sys.argv[1]
dir_list = ["/users/asebe/bvanwaardenburg/git/data/241204_ParseBio_Nvec_Tcas_nuclei/Tcas/mapping_STARsolo/mapping_STARsolo_N/", 
            "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas_sep/mapping_STARsolo/mapping_STARsolo_N/"]

base_dir = "/users/asebe/bvanwaardenburg/git/data/"

file_list = []
pattern = "**/Solo.out/Gene/raw/matrix.mtx"
for paths in dir_list:
    all_files = glob.glob(os.path.join(paths, pattern), recursive=True)
    file_list.append(all_files)
file_list.sort(reverse=True)

# Extracts the full basename
def extract_basename(file_path):
    """
    Returns basename:
      'BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_DSP'
    from a path:
      '/.../BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_DSP/Solo.out/Gene/raw/matrix.mtx'
    """
    return os.path.basename(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(
                    os.path.dirname(file_path)
                    )
            )
        )
    )

def simplify_name(base_name):
    label = base_name.split("_")[0]         # e.g. 'BCA001' from 'BCA001_lib_..._DSP'
    parts = base_name.split("_")            # e.g. 'DSP' from 'BCA001_lib_..._DSP'
    return label + "_" + "_".join(parts[4:])


# =================  Knee Plot ================= #

# ------------------------------------------------------------------
# PLOT SETTINGS
# ------------------------------------------------------------------
# Create a color palette—one distinct color per file
# colors = sns.color_palette(["#B3B3B3", "#C1B7AA", "#CEBCA1", "#DCC198", "#E8C687", "#EFCC6B", "#FDD833", "#ABD15C", "#BDBB7B", "#E190BA", "#919FC6",
#                            "#B099A9", "#CF948C", "#EE8F6E", "#E4956C", "#BAA47F", "#90B392", "#66C2A5"]) 
colors = sns.color_palette(["#F6D24F", "#E9D837", "#D0D842", "#B7D84C", "#CFA69B", "#D58EC4", "#BB94C6", "#A29AC9",
                            "#E4956C", "#BAA47F", "#90B392", "#66C2A5"]) 
colors += colors

# ------------------------------------------------------------------
# PLOT: Combined Knee Plot
# ------------------------------------------------------------------
plt.figure(figsize=(10, 6))
count = 0

for dir_path in file_list:
    dir_path.sort(reverse=True)
    for file_path in dir_path:

        sample_label = extract_basename(file_path)  # e.g. 'BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_DSP'
        sample_label = simplify_name(sample_label)

        # Convert to sparse CSC format
        sparse_matrix = scipy.io.mmread(file_path).tocsc()  

        # Sum counts per barcode (column-wise sum) and sort
        barcode_counts = np.array(sparse_matrix.sum(axis=0)).flatten()
        barcode_counts_sorted = np.sort(barcode_counts)[::-1]

        # Generate rank (x-axis)
        ranks = np.arange(1, len(barcode_counts_sorted) + 1)

        plt.plot(
            ranks,
            barcode_counts_sorted,
            label=sample_label,
            linewidth=1,
            color=colors[count]
        )
        count += 1

plt.xscale("log")
plt.yscale("log")
plt.xlabel("# Barcodes (log scale)")
plt.ylabel("# Transcripts (log scale)")
plt.title("Combined knee plot across Tcas samples")
plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.tight_layout()
plt.savefig("%s/knee_plots_tcas.png" %(base_dir), dpi=300)
plt.close()

