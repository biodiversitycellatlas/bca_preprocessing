#!/usr/bin/env python

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
souporcell_dir = "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF/souporcell/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM_12clusters/"
pca_df = pd.read_csv(f"{souporcell_dir}clusters.tsv", sep="\t")

# Modify the assignment column
pca_df['assignment'] = pca_df.apply(
    lambda row: "NA" if row['status'] == "unassigned" or '/' in row['assignment'] else row['assignment'],
    axis=1
)

# Select columns starting with 'c' and assignment
cols_to_keep = [col for col in pca_df.columns if col.startswith('c')]
cols_to_keep.append('assignment')
pca_df = pca_df[cols_to_keep]

# Max normalize the first 4 columns
pca_df.iloc[:, :4] = pca_df.iloc[:, :4].div(pca_df.iloc[:, :4].max(axis=1), axis=0)

# Perform PCA
features = pca_df.iloc[:, :4] 
pca = PCA(n_components=min(features.shape[0], features.shape[1]))
pca_result = pca.fit_transform(features)

# Create a df for PCA results
res_pca_df = pd.DataFrame(pca_result[:, :2], columns=['Dim.1', 'Dim.2'])
res_pca_df['assignment'] = pca_df['assignment']

# PCA plot
plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=res_pca_df,
    x='Dim.1',
    y='Dim.2',
    hue='assignment',
    s=10
)
plt.title("PCA Plot of Cells by Genotype Scores")
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.legend(title="Assignment", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("PCA_paper.png")

