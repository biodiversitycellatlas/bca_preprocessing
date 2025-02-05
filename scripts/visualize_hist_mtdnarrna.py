#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import glob


# ------------------------------------------------------------------
# LOAD FILES
# ------------------------------------------------------------------
filepaths = ["/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002/mapping_stats.tsv", 
                  "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA003_BCA004_FINAL/mapping_stats.tsv",
                  "/users/asebe/bvanwaardenburg/git/data/241204_ParseBio_Nvec_Tcas_nuclei/Nvec_FINAL/mapping_stats.tsv"]
# dir_list = sys.argv[1]

base_dir = "/users/asebe/bvanwaardenburg/git/data/"

# Columns of interest
columns_of_interest = [
    "% rRNA in Unique reads",
    "%rRNA in multimappers all pos",
    "%rRNA in multimappers primary pos",
    "% mtDNA in Unique reads",
    "%mtDNA in multimappers all pos",
    "%mtDNA in multimappers primary pos"
]

# Dictionary to store data per column across all samples
data_dict = {col: [] for col in columns_of_interest}
sample_labels = {col: [] for col in columns_of_interest}

# Read data from each TSV file
for filepath in filepaths:
    try:
        df = pd.read_csv(filepath, sep='\t')

        # Filter rows where 'Directory' column ends with '_N'
        df_filtered = df[df['Directory'].str.endswith('_N', na=False)]

        if not df_filtered.empty:
            for _, row in df_filtered.iterrows():
                base_name = row['Sample']
                label = base_name.split("_")[0]  # e.g. 'BCA001' from 'BCA001_lib_..._DSP'
                parts = base_name.split("_")[4:]  # e.g. 'DSP' from 'BCA001_lib_..._DSP'
                cond = "_".join(parts)
                simplified_label = f"{label}_{cond}"

                for col in columns_of_interest:
                    if col in df_filtered.columns:
                        data_dict[col].append(row[col]*100)
                        sample_labels[col].append(simplified_label) 

    except Exception as e:
        print(f"Error processing {filepath}: {e}")


# Assign consistent colors to labels
unique_labels = list(set(label for labels in sample_labels.values() for label in labels))
color_map = {label: plt.cm.tab20(i % len(unique_labels)) for i, label in enumerate(unique_labels)}  # Use tab10 colormap for distinct colors

# Create a single figure with subplots
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(18, 12))
axes = axes.flatten()

for i, col in enumerate(columns_of_interest):
    ax = axes[i]
    colors = [color_map[label] for label in sample_labels[col]]

    ax.bar(sample_labels[col], data_dict[col], color=colors)
    ax.set_xlabel("Sample")
    ax.set_ylabel(col)
    ax.set_title(f"{col} ")
    ax.tick_params(axis='x', rotation=90)

plt.tight_layout()
plt.savefig("%s/hist_plots_mm.png" %(base_dir), dpi=300)
plt.close()
