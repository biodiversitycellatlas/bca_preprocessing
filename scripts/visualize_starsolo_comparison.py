#!/usr/bin/env python

import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


# ------------------------------------------------------------------
# LOAD FILES
# ------------------------------------------------------------------
# directory_path = "/users/asebe/bvanwaardenburg/git/data/241106_BD_Rhapsody_Nvec/Nvec_FINAL/mapping_STARsolo/mapping_STARsolo_CR/"
directory_path = sys.argv[1]

pattern = "**/Solo.out/Gene/UMIperCellSorted.txt"
all_files = glob.glob(os.path.join(directory_path, pattern), recursive=True)
all_files.sort()

# Extracts the full basename
def extract_basename(file_path):
    """
    Returns basename:
      'BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_DSP'
    from a path:
      '/.../BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_DSP/Solo.out/Gene/UMIperCellSorted.txt'
    """
    return os.path.basename(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(file_path)
            )
        )
    )


# ------------------------------------------------------------------
# BUILD A DICTIONARY
# Will create a nested dictionary with the following structure:
#   files_dict[label][condition] = file_path
#   files_dict[label][condition] = file_path
# ------------------------------------------------------------------
files_dict = {}
labels = []
no_conditions = False

# Populate files_dict
for fpath in all_files:
    bname = extract_basename(fpath)     # e.g. 'BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_DSP'
    label = bname.split("_")[0]         # e.g. 'BCA001' from 'BCA001_lib_..._DSP'
    parts = bname.split("_")            # e.g. 'DSP' from 'BCA001_lib_..._DSP'
    condition = "_".join(parts[4:])
    print(parts)

    if label in files_dict:
        files_dict[label][condition] = fpath
    elif len(parts) < 4:
        files_dict[label] = fpath
        no_conditions = True
    elif label not in files_dict:
        labels.append(label)
        files_dict[label] = {} 
        files_dict[label][condition] = fpath


# ------------------------------------------------------------------
# DETERMINE COMMON CONDITIONS
# ------------------------------------------------------------------
if no_conditions == False:
    bca1_conditions = set(files_dict[labels[0]].keys())
    bca2_conditions = set(files_dict[labels[1]].keys())
    print(bca1_conditions, bca2_conditions)

    # Compute common and unique conditions
    common_conditions = sorted(bca1_conditions & bca2_conditions)  # Intersection
    unique_conditions = sorted(bca1_conditions ^ bca2_conditions)  # Symmetric difference
    all_conditions = common_conditions + unique_conditions  # Combine them
else:
    all_conditions = set(files_dict.keys())

print("All conditions:", all_conditions)


# ================= 1) Overlaid Knee Plot ================= #

# ------------------------------------------------------------------
# PLOT SETTINGS
# ------------------------------------------------------------------
# Create a color palette—one distinct color per file
colors = sns.color_palette("tab10", len(all_conditions)) 
colors += colors

# ------------------------------------------------------------------
# PLOT: Overlaid Knee Plot
# ------------------------------------------------------------------
plt.figure(figsize=(8, 6))

for i, file_path in enumerate(all_files):
    
    sample_label = extract_basename(file_path)  # e.g. 'BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_DSP'
    umi_data = pd.read_csv(file_path, sep="\t", header=None, names=["UMICount"])
    
    # Sort descending & Rank values
    umi_data.sort_values("UMICount", ascending=False, inplace=True, ignore_index=True)
    umi_data["Rank"] = np.arange(1, len(umi_data) + 1)

    plt.plot(
        umi_data["Rank"],
        umi_data["UMICount"],
        label=sample_label,
        color=colors[i],
        linewidth=1
    )

plt.xscale("log")
plt.yscale("log")
plt.xlabel("# Barcodes (log scale)")
plt.ylabel("# Transcripts (log scale)")
plt.title("Combined knee plot")
plt.legend(loc="best")

plt.tight_layout()
plt.savefig("%s/knee_plot_overlay.png" %(directory_path), dpi=300)
plt.close()



# ================= 2) Distribution Plots Stacked ================= #

# ------------------------------------------------------------------
# PLOT SETTINGS
# ------------------------------------------------------------------
threshold = 10
color_1 = "#0077b6" 
color_2 = "#e76f51"  

fig, axes = plt.subplots(nrows=len(all_conditions), ncols=1,
                         figsize=(10, 3 * len(all_conditions)))

# If there's only one condition, axes is not a list
if len(all_conditions) == 1:
    axes = [axes]

# ------------------------------------------------------------------
# PLOT: Distribution Plots Stacked
# ------------------------------------------------------------------
for i, cond in enumerate(all_conditions):
    ax = axes[i]
    if no_conditions == False:
        # Retrieve BCA001 and BCA002 file paths
        file_path_1 = files_dict[labels[0]][cond]
        file_path_2 = files_dict[labels[1]][cond]
        
        # Extract sample labels from the path
        sample_label_1 = extract_basename(file_path_1)
        sample_label_2 = extract_basename(file_path_2)
        
        # ---- Read data
        umi_data_1 = pd.read_csv(file_path_1, sep="\t", header=None, names=["UMICount"])
        umi_data_2 = pd.read_csv(file_path_2, sep="\t", header=None, names=["UMICount"])

        # ---- Filter data
        filtered_1 = umi_data_1[umi_data_1["UMICount"] > threshold]
        filtered_2 = umi_data_2[umi_data_2["UMICount"] > threshold]

        # ---- Plot both
        sns.histplot(
            filtered_1["UMICount"],
            bins=50, kde=True,
            ax=ax,
            log_scale=(True, False),
            color=color_1,
            alpha=0.4,
            label=sample_label_1
        )
        sns.histplot(
            filtered_2["UMICount"],
            bins=50, kde=True,
            ax=ax,
            log_scale=(True, False),
            color=color_2,
            alpha=0.4,
            label=sample_label_2
        )
    else:
        # Retrieve BCA001 and BCA002 file paths
        file_path = files_dict[cond]
        
        # Extract sample labels from the path
        sample_label = extract_basename(file_path)
        
        # ---- Read data
        umi_data = pd.read_csv(file_path, sep="\t", header=None, names=["UMICount"])

        # ---- Filter data
        filtered = umi_data[umi_data["UMICount"] > threshold]

        # ---- Plot both
        sns.histplot(
            filtered["UMICount"],
            bins=50, kde=True,
            ax=ax,
            log_scale=(True, False),
            color=color_1,
            alpha=0.4,
            label=sample_label
        )
    
    ax.legend()
    ax.set_title(f"{cond}")
    ax.set_xlabel("UMI / Cell Barcode (Log Scale, UMI > %s)" %(threshold))
    ax.set_ylabel("Density")

plt.tight_layout()
plt.savefig(f"{directory_path}/distribution_plots_stacked.png", dpi=300)
plt.close()
