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
# directory_path = sys.argv[1]
directory_path = "/users/asebe/bvanwaardenburg/git/data/250115_ParseBio_Nvec_Tcas_Pliv_Cele/Nvec_BCA009_BCA010/mapping_STARsolo/mapping_STARsolo_CRNF_SPIPE/"

pattern = "**/*_Solo.out/Gene/UMIperCellSorted.txt"
all_files = glob.glob(os.path.join(directory_path, pattern), recursive=True)
all_files.sort()

# Extracts the full basename
def extract_basename(file_path):
    """
    Returns basename, e.g.:
      'BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACME_GM'
    from a path:
      '/.../BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACME_GM/Solo.out/Gene/UMIperCellSorted.txt'
    """
    return os.path.basename(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(file_path)
            )
        )
    )

files_dict = {}
labels = []         # Now this will be a list of fixation techniques
no_conditions = False

# ------------------------------------------------------------------
# BUILD A DICTIONARY
#   files_dict[fixation_technique][library] = file_path
# ------------------------------------------------------------------
for fpath in all_files:
    bname = extract_basename(fpath)
    parts = bname.split("_")
    # Example: ["BCA001", "lib", "13077AAF", "CAGATCAC-ATGTGAAG", "ACME", "GM"]

    # 1) "label" = fixation technique (last part)
    # 2) "condition" = library ID (the first part)
    label = "_".join(parts[4:])     # e.g. "ACME_GM"
    condition = parts[0]           # e.g. "BCA001"

    if label not in files_dict:
        labels.append(label)
        files_dict[label] = {}

    # If for some reason the name has fewer than 4 parts, handle no_conditions:
    if len(parts) < 5:
        files_dict[label] = fpath
        no_conditions = True
    else:
        files_dict[label][condition] = fpath


# ------------------------------------------------------------------
# DETERMINE COMMON CONDITIONS (now libraries)
# ------------------------------------------------------------------
if not no_conditions:
    # For simplicity, assume exactly two libraries so the old logic works
    # (BCA001, BCA002).  If you have more libraries, this part needs
    # generalizing—currently it only handles the first two.
    libset1 = set(files_dict[labels[0]].keys())
    libset2 = set(files_dict[labels[1]].keys())

    common_conditions = sorted(libset1 & libset2)  # intersection
    unique_conditions = sorted(libset1 ^ libset2)  # symmetric difference
    all_conditions = common_conditions + unique_conditions
else:
    # If we had no conditions, we just store the keys
    all_conditions = set(files_dict.keys())

print("All conditions (now referring to libraries):", all_conditions)
print("All labels (now referring to fixation techniques):", labels)

# ================= 1) Overlaid Knee Plot ================= #

# Create a color palette—one distinct color per file
colors = sns.color_palette("tab10", len(all_files))
colors += colors

plt.figure(figsize=(8, 6))

for i, file_path in enumerate(all_files):
    sample_label = extract_basename(file_path)
    umi_data = pd.read_csv(file_path, sep="\t", header=None, names=["UMICount"])
    
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
plt.savefig(f"{directory_path}/knee_plot_overlay.png", dpi=300)
plt.close()


# ================= 2) Distribution Plots Stacked ================= #

threshold = 10
color_1 = "#0077b6" 
color_2 = "#e76f51"

fig, axes = plt.subplots(nrows=len(all_conditions), ncols=1,
                         figsize=(10, 3 * len(all_conditions)))

# If there's only one condition, axes is not a list
if len(all_conditions) == 1:
    axes = [axes]

for i, cond in enumerate(all_conditions):
    ax = axes[i]
    if not no_conditions:
        # Retrieve the two fixation techniques from 'labels' 
        # (which we are now using for ACME_GM, DSP, etc.)
        file_path_1 = files_dict[labels[0]].get(cond, None)
        file_path_2 = files_dict[labels[1]].get(cond, None)
        
        # If either is None, it means that library doesn't exist for that technique.
        if file_path_1:
            umi_data_1 = pd.read_csv(file_path_1, sep="\t", header=None, names=["UMICount"])
            filtered_1 = umi_data_1[umi_data_1["UMICount"] > threshold]
            sns.histplot(
                filtered_1["UMICount"],
                bins=50, kde=True,
                ax=ax,
                log_scale=(True, False),
                color=color_1,
                alpha=0.4,
                label=f"{labels[0]}_{cond}"
            )
        if file_path_2:
            umi_data_2 = pd.read_csv(file_path_2, sep="\t", header=None, names=["UMICount"])
            filtered_2 = umi_data_2[umi_data_2["UMICount"] > threshold]
            sns.histplot(
                filtered_2["UMICount"],
                bins=50, kde=True,
                ax=ax,
                log_scale=(True, False),
                color=color_2,
                alpha=0.4,
                label=f"{labels[1]}_{cond}"
            )
    else:
        # No condition split
        file_path = files_dict[cond]
        umi_data = pd.read_csv(file_path, sep="\t", header=None, names=["UMICount"])
        filtered = umi_data[umi_data["UMICount"] > threshold]
        sns.histplot(
            filtered["UMICount"],
            bins=50, kde=True,
            ax=ax,
            log_scale=(True, False),
            color=color_1,
            alpha=0.4,
            label=cond
        )
    
    ax.legend()
    ax.set_title(f"Library: {cond}")
    ax.set_xlabel(f"UMI / Cell Barcode (Log Scale, UMI > {threshold})")
    ax.set_ylabel("Density")

plt.tight_layout()
plt.savefig(f"{directory_path}/distribution_plots_stacked.png", dpi=300)
plt.close()
