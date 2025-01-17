import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
file_path = "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF/mapping_STARsolo_N_addattrib/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM/Solo.out/Gene/UMIperCellSorted.txt"
umi_data = pd.read_csv(file_path, sep="\t", header=None, names=["UMICount"])

# Generate rank for the x-axis
umi_data["Rank"] = np.arange(1, len(umi_data) + 1)

# Define tscp_cutoff and compute statistics
tscp_cutoff = 294
median_tscp = umi_data["UMICount"].median()
num_cells = (umi_data["UMICount"] >= tscp_cutoff).sum()

# Color points based on tscp_cutoff
umi_data["Color"] = np.where(umi_data["UMICount"] >= tscp_cutoff, "green", "gray")

# Create the plot
plt.figure(figsize=(8, 6))
plt.scatter(umi_data["Rank"], umi_data["UMICount"], c=umi_data["Color"], alpha=0.6, s=10)

plt.xscale("log")
plt.yscale("log")
plt.xlabel("# Barcodes (logscale)")
plt.ylabel("# Transcripts (logscale)")
plt.title("Identified Cells")

# Add annotations in the bottom-left corner 
text = (
    f"Number of Cells: {num_cells}\n"
    f"Transcript Cutoff: {tscp_cutoff}\n"
    f"Median Transcripts: {median_tscp:.0f}"
)
plt.text(
    0.05, 0.05,
    text,
    transform=plt.gca().transAxes,
    fontsize=10,
    color="black"
)


plt.savefig("knee_plot.png")
plt.close()
