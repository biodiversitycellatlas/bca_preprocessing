#!/usr/bin/env python

import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import re


# ------------------------------------------------------------------
# LOAD FILES
# ------------------------------------------------------------------
# file_path = "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002/mapping_STARsolo/intron_vs_mtdna_metrics.csv"
# file_path = sys.argv[1]

directory_path = "/users/asebe/bvanwaardenburg/git/data/"
pattern = "**/mapping_STARsolo/intron_vs_mtdna_metrics.csv"
all_files = glob.glob(os.path.join(directory_path, pattern), recursive=True)
all_files.sort()

colors = sns.color_palette("Spectral", 30) 
c = 0

plt.figure(figsize=(12, 8))
for n, file_path in enumerate(all_files):
    
    df = pd.read_csv(file_path)   
    df[["label", "condition"]] = df["sample_name"].str.extract(r"(BCA00\d+)_.*?_[ACTG]{8}-[ACTG]{8}_(.+)$")
    df["sample_name"] = df["label"] + "_" + df["condition"]

    # Convert numeric columns to float
    df['mtDNA_fraction'] = pd.to_numeric(df['mtDNA_fraction'])
    df["mtDNA_fraction"] *= 100
    df['fraction_intronic'] = pd.to_numeric(df['fraction_intronic'])
    df["fraction_intronic"] *= 100

    # Filter for 'N' config & set negative values to zero
    df["fraction_intronic"] = df["fraction_intronic"].clip(lower=0)
    df_filtered = df[df["config"] == "N"]
    
    # Plot
    for i, sample in enumerate(df_filtered["sample_name"]):
        sample_data = df_filtered[df_filtered["sample_name"] == sample]
        plt.scatter(sample_data["mtDNA_fraction"], sample_data["fraction_intronic"], label=sample_data["sample_name"].iloc[0], alpha=0.7, color=colors[c])
        c +=1

plt.xlabel("Percentage of mtDNA in Unique reads")
plt.ylabel("Fraction of UMIs mapped to introns")
plt.title("Comparing the percentage of mtDNA versus Intronic counts across samples")
plt.legend(title="Sample Name", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.grid(True)
plt.tight_layout()
plt.savefig("mtDNA_vs_intron.png", dpi=300)
plt.close()

