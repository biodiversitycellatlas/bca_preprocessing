#!/usr/bin/env Rscript
# =============================================================================
# This script searches the base directory recursively for mapping results,
# processes the STARsolo output, and creates visualizations.
# =============================================================================

# =============================================================================
# Load Packages
# =============================================================================

options(repos = c(CRAN = "https://cran.r-project.org"))
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('Matrix')) install.packages('Matrix'); library('Matrix')
if (!require('ggridges')) install.packages('ggridges'); library('ggridges')
if (!require('scales')) install.packages('scales'); library('scales')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('patchwork')) install.packages('patchwork'); library('patchwork')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')


# =============================================================================
# Load arguments and set up directories
# =============================================================================

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript create_visualizations.R <data_directory>")
}
data_dir <- args[1]

# Define a base font size and set a global theme 
base_font <- 14
theme_set(
  theme_minimal(base_size = base_font) +
    theme(
      axis.title     = element_text(size = base_font),
      axis.text      = element_text(size = base_font * 0.9),
      plot.title     = element_text(size = base_font * 1.2, face = "bold"),
      strip.text     = element_text(size = base_font),
      legend.text    = element_text(size = base_font),
      legend.title   = element_text(size = base_font)
    )
)

# Path to STARsolo outputs 
mapping_dir <- file.path(data_dir, "mapping_STARsolo")

# Find all sample‐ids by listing first‐level dirs
sample_ids <- list.dirs(mapping_dir, recursive = FALSE, full.names = FALSE)
raw_dirs   <- list()

for (sample_id in sample_ids) {
  path_pattern   <- file.path(mapping_dir, sample_id, "*_Solo.out", "GeneFull_Ex50pAS", "raw")
  expanded_paths <- Sys.glob(path_pattern)
  if (length(expanded_paths) > 0) {
    # store vector of raw paths under sample_id
    raw_dirs[[sample_id]] <- expanded_paths
  }
}

# =============================================================================
# Function to read STARsolo data
# =============================================================================
read_star_data <- function(directory, dataset_name) {
  # Read matrix, features, and barcodes
  mat <- Matrix::readMM(file.path(directory, "matrix.mtx"))
  features <- data.table::fread(file.path(directory, "features.tsv"), header = FALSE)
  barcodes_dt <- data.table::fread(file.path(directory, "barcodes.tsv"), header = FALSE)
  
  # Assign row and column names to the matrix
  rownames(mat) <- features$V1
  colnames(mat) <- barcodes_dt$V1
  
  # Compute per-cell summaries
  cell_umis <- Matrix::colSums(mat)
  cell_genes <- Matrix::colSums(mat > 0)
  
  # Create a data.table for cells
  cell_dt <- data.table(
    cells      = colnames(mat),
    cell_sizes = cell_umis,
    cell_genes = cell_genes
  )
  
  # Compute per-gene summaries
  gene_counts <- Matrix::rowSums(mat)
  gene_dt <- data.table(
    gene      = rownames(mat),
    gene_umis = gene_counts
  )
  
  return(list(cell_summary = cell_dt, gene_summary = gene_dt))
}

# =============================================================================
# Process each dataset
# =============================================================================
cell_summary_list <- list()
gene_summary_list <- list()

# For each raw‐folder under this sample_id and store with simple sample_id key
for (dataset_name in names(raw_dirs)) {
  for (directory in raw_dirs[[dataset_name]]) {
    res <- read_star_data(directory, dataset_name)
    cell_summary_list[[dataset_name]] <- res$cell_summary
    gene_summary_list[[dataset_name]] <- res$gene_summary
  }
}

# =============================================================================
# Combine & Filter Data
# =============================================================================
cell_summary_dt <- rbindlist(cell_summary_list, idcol = "dataset")
gene_summary_dt <- rbindlist(gene_summary_list, idcol = "dataset")

# Define thresholds for filtering
umi_thr <- c(10, 10000)  # UMI count thresholds per cell
gen_thr <- 10           # Minimum genes detected per cell

# Filter cells based on UMI and gene count thresholds
filt_cells <- cell_summary_dt[
  cell_sizes >= umi_thr[1] & cell_sizes <= umi_thr[2] & cell_genes >= gen_thr
]

# Check if any cells pass the filters
if (nrow(filt_cells) == 0)
    stop("No cells pass the filters, adjust umi_thr or gen_thr")

# Filter cells and genes for plotting (only positive counts)
sum_cells_plot <- cell_summary_dt[cell_sizes > 0]
sum_genes_plot <- gene_summary_dt[gene_umis > 0]

# Define a dynamic color palette
dataset_names <- unique(cell_summary_dt$dataset)
num_datasets <- length(dataset_names)
data_cols <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(num_datasets), dataset_names)

# Remove duplicates if any
sum_cells_plot <- sum_cells_plot %>% distinct()
filt_cells <- filt_cells %>% distinct()

# =============================================================================
# UMI Distribution Plot (Top 10,000 cells per dataset)
# =============================================================================
cell_n <- 10000
ttl <- sprintf("Top %s BC", cell_n)
sub_label <- sprintf("top%s", cell_n)

dist_dt <- sum_cells_plot[order(dataset, -cell_sizes), .SD[1:pmin(.N, cell_n)], by = dataset]
dist_dt <- dist_dt[cell_sizes > 0]

gp_umis_dist <- ggplot(dist_dt, aes(cell_sizes, dataset, fill = dataset)) +
  ggridges::geom_density_ridges2(scale = 0.9, rel_min_height = 0.01) +
  geom_vline(xintercept = 345, linetype = "dashed", color = "red") +
  scale_x_continuous(
    breaks = 10^seq(1, 5),
    trans = "log10", labels = scales::trans_format('log10', scales::math_format(10^.x))
  ) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = data_cols) +
  theme(
    strip.text = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "", x = "UMIs / cell barcode", title = ttl)

gp_cells <- ggplot(filt_cells, aes(dataset, fill = dataset)) +
  geom_bar(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = data_cols) +
  geom_text(aes(label = ..count..), stat = "count", vjust = 0.5, hjust = -0.2) +
  coord_flip() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = "", title = sprintf("Cells\n(%s<UMI<%s)", umi_thr[1], umi_thr[2]))

# Create output directory for figures
scdb_fig_dir <- "R_images"
map_dir <- basename(data_dir)  # use the base directory name for labeling

if (!dir.exists(scdb_fig_dir)) {
  dir.create(scdb_fig_dir, recursive = TRUE)
}

# Combine plots with relative widths
combined_plot <- gp_umis_dist + gp_cells +
  plot_layout(nrow = 1)

# Save the combined_plot
ggsave(file.path(scdb_fig_dir, sprintf("UMI_dist_%s_%s.png", map_dir, sub_label)), 
       combined_plot, dpi=300, width  = 14, height = 5 )

# =============================================================================
# Per Dataset & Per Cell Plots
# =============================================================================
gp_cells <- ggplot(filt_cells, aes(dataset, fill = dataset)) +
  geom_bar(color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = data_cols) +
  geom_text(aes(label = ..count..), stat = "count", vjust = 0.5, angle = 90, hjust = -0.2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  labs(x = "", y = "", title = sprintf("Cells\n(%s<UMI<%s)", umi_thr[1], umi_thr[2]))

gp_genes <- ggplot(unique(sum_genes_plot[gene_umis > 0, .(dataset, gene)]), aes(dataset, fill = dataset)) +
  geom_bar(color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = data_cols) +
  geom_text(aes(label = ..count..), stat = "count", vjust = 0.5, angle = 90, hjust = -0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") +
  labs(x = "", y = "", title = "Genes\n(UMI>0)")

gp_umis_cells <- ggplot(filt_cells, aes(dataset, cell_sizes, fill = dataset)) +
  geom_boxplot() +
  geom_text(data = filt_cells[, median(cell_sizes), .(dataset)], aes(label = round(V1), y = V1), hjust = 0.5, vjust = -0.5) +
  scale_y_continuous(trans = "log10", limits = c(umi_thr[1], NA)) +
  scale_fill_manual(values = data_cols) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  labs(x = "", y = "", title = "UMIs / cell")

gp_genes_cells <- ggplot(filt_cells, aes(dataset, cell_genes, fill = dataset)) +
  geom_boxplot() +
  geom_text(data = filt_cells[, median(cell_genes), .(dataset)], aes(label = round(V1), y = V1), hjust = 0.5, vjust = -0.5) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = data_cols) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  labs(x = "", y = "", title = "Genes / cell")

gp_umi_genes <- ggplot(sum_genes_plot, aes(dataset, gene_umis, fill = dataset)) +
  geom_boxplot() +
  geom_text(data = sum_genes_plot[gene_umis > 0, median(gene_umis), .(dataset)], aes(label = round(V1), y = V1), hjust = 0.5, vjust = -0.5) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") +
  scale_fill_manual(values = data_cols) +
  labs(x = "", y = "", title = "UMIs / gene")

gps <- (gp_cells / gp_genes) | (gp_umis_cells / gp_genes_cells / gp_umi_genes)
ggsave(file.path(scdb_fig_dir, sprintf("cells_genes_%s.png", map_dir)), gps, width = 26, height = 15)
