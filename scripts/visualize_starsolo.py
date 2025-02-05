import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Load data
output_dir = "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002/mapping_STARsolo/mapping_STARsolo_N/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM"
file_path = output_dir + "/Solo.out/Gene/UMIperCellSorted.txt"
umi_data = pd.read_csv(file_path, sep="\t", header=None, names=["UMICount"])

# Create rank for the x-axis
umi_data["Rank"] = np.arange(1, len(umi_data) + 1)

# Define tscp_cutoff and compute statistics
tscp_cutoff = 294
# 4. Compute statistics
median_tscp = umi_data["UMICount"].median()
num_cells = (umi_data["UMICount"] >= tscp_cutoff).sum()

# 5. Assign colors based on the cutoff
umi_data["Color"] = np.where(umi_data["UMICount"] >= tscp_cutoff, "green", "gray")


# Create scatterplot
plt.figure(figsize=(8, 6))
plt.scatter(umi_data["Rank"], umi_data["UMICount"], c=umi_data["Color"], alpha=0.6, s=10)
plt.xlim(0, 8000)
#plt.xscale("log")
plt.yscale("log")
plt.ylim( (pow(10,1),pow(10,4)) )
plt.xlabel("Sorted Barcode ID")
plt.ylabel("Unique UMI counts")
plt.title("UMI distribution plot for BCA001_ACMEsorb_GM")
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
plt.savefig("%s/knee_plot.png" %(output_dir), dpi=300)
plt.close()


# Second plot 
threshold = 50
filtered_data = umi_data[umi_data["UMICount"] > threshold]

plt.figure(figsize=(10, 6))
sns.histplot(filtered_data["UMICount"], bins=50, kde=True, log_scale=(True, False))
plt.title("Distribution of Read Counts per Cell Barcode (Filtered)")
plt.xlabel("Read Counts per CB (Log Scale)")
plt.ylabel("Density")
plt.savefig("%s/knee_plot_2.png" %(output_dir), dpi=300)
plt.close()
