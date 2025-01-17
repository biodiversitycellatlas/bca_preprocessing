#!/usr/bin/env python

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import upsetplot

# Import data and load it into a df
base_dir = "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF/mapping_STARsolo_N_addattrib"

single_file = base_dir + "/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM/sorted_multimapper_gene_counts.txt"
single_data = pd.read_csv(single_file, sep="\t", header=None, names=["GeneID", "Count"])

# Loop through all directories in the base_dir and append the file names to a list
comparison_files = []
for base, dirs, files in os.walk(base_dir):
    for file in files:
        if file == "sorted_multimapper_gene_counts.txt":
            comparison_files.append(os.path.join(base, file))

print(comparison_files)

comparison_data = []
for file in comparison_files:
    # Get the base name of the file
    base_name = os.path.basename(os.path.dirname(file))

    # Read the file contents
    df = pd.read_csv(file, sep="\t", header=None, names=["GeneID", "Count"])

    # Append the base name and the columns to comparison_data
    for _, row in df.iterrows():
        comparison_data.append([base_name, row["GeneID"], row["Count"]])

# Convert to DataFrame for easier handling
comparison_df = pd.DataFrame(comparison_data, columns=["BaseName", "GeneID", "Count"])

# Display the DataFrame
print(comparison_df.head())


# Remove zero counts to avoid the warning during variance calculation
filtered_data = single_data[single_data["Count"] > 0]  

# Histogram of log scaled number of counts across genes
plt.figure(figsize=(10, 6))
sns.histplot(filtered_data["Count"], bins=30, log_scale=True, kde=True)
plt.xlabel("Read Counts (log scale)")
plt.ylabel("Frequency")
plt.title("Diversity of gene counts among multi-mapped reads (log scaled)")
plt.savefig("count_histogram.png", dpi=300)
plt.close()


# Select the top 20 genes with the highest counts
top_20_genes = filtered_data.nlargest(20, "Count")

# Create a bar plot for the top 20 genes with the highest counts
sns.barplot(y=top_20_genes["GeneID"], x=top_20_genes["Count"])
plt.xlabel("Read Counts")
plt.ylabel("Gene ID")
plt.title("Top 20 genes with the highest gene counts among multi-mapped reads")
plt.tight_layout()
plt.savefig("top_20_genes_histogram.png", dpi=300)
plt.close()


# Count occurrences of each GeneID in BaseName
gene_base_counts = comparison_df.groupby(["GeneID", "BaseName"])["Count"].sum().reset_index(name="Occurrences")

# Identify the top 20 genes by total count
top_genes = comparison_df.groupby("GeneID")["Count"].sum().nlargest(20).index

# Filter the data for only the top 20 genes
filtered_gene_base_counts = gene_base_counts[gene_base_counts["GeneID"].isin(top_genes)]

# Create a pivot table for the filtered data
filtered_heatmap_data = filtered_gene_base_counts.pivot(
    index="GeneID", columns="BaseName", values="Occurrences"
).fillna(0)

# Plot heatmap for the top 20 genes
plt.figure(figsize=(12, 8))
sns.heatmap(
    filtered_heatmap_data,
    annot=True,
    fmt="g",
    cmap="coolwarm",
    cbar_kws={'label': 'Counts'},
    linewidths=0.5,
    linecolor='gray'
)
plt.title("Top 20 most frequent genes among different samples in multi-mapped reads", fontsize=14)
plt.xlabel("Sample", fontsize=10)
plt.ylabel("GeneID", fontsize=10)
plt.xticks(fontsize=8, rotation=45, ha='right') 
plt.yticks(fontsize=8)
plt.tight_layout()

plt.savefig("heatmap_gene_occurences.png", dpi=300)
plt.close()

