# libraries
library("metacell")
source("/ext_programs/metacell-downstream-functions/scripts/Downstream_functions.R")
source("/ext_programs/metacell-downstream-functions/scripts/helper.R")
source("/ext_programs/metacell-downstream-functions/scripts/Modified_functions.R")

# metacell parameters
tgconfig::set_param("mc_cores", 1, "metacell") # this allows gstat to work without crashing due to excessive memory usage (and it's equally fast)

# run name (objects in the metacell db will be given this suffix)
run_name = "nvec"

# init metacell database
# at the beginning, there's just a `mat` object
metacell::scdb_init("/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002/\
    mapping_STARsolo_N/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM/Solo.out/GeneFull/filtered",force_reinit=TRUE)


#### Load data ####

# quality control figures wil go here
metacell::scfigs_init("/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Metacell/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM/")

# load umi table
mat = metacell::scdb_mat("nvec")

# cell sizes
mat_cz = Matrix::colSums(mat@mat)



#### Cell filtering ####

# select cells that are too small or too large
# thresholds
cz_large_thr = 15000
cz_small_thr = 100
# apply them
cells_large = names(which(mat_cz > cz_large_thr))
cells_small = names(which(mat_cz < cz_small_thr))

# keep cells that are of intermediate size
cells_keep = names( which(mat_cz >= cz_small_thr & mat_cz <= cz_large_thr) )

# visualise cell size distribution to see if the above thresholds make sense
pdf(sprintf("results_figs/%s.cs_distribution.pdf", "nvec"), height = 6, width = 5)
# plot histogram
hh = hist(
	log10(mat_cz + 1),
	breaks = 60,
	main = sprintf("UMI/cell %s", "nvec"),
	col = "gray", 
	border = NA,
	xlab = "log10 UMI+1",
	xlim = c(0,5),
	las = 1)
hist(log10(mat_cz + 1) [ names(mat_cz) %in% cells_keep ], col = alpha("blue", 0.7), add = TRUE, breaks = 60, border = NA)
legend(
	"topleft", 
	legend = c(
		sprintf("all cells n = %i", length(mat_cz)), 
		sprintf("kept cells n = %i", length(cells_keep))), 
	fill = c("gray", "blue", "magenta3"), bty = "n", cex = 0.7, border = NA)
# title(sub = sprintf("n = %i / %i / %i cells in each interval (out of %i)", length(cells_small), length(cells_keep), length(cells_large), length(mat_cz)), cex.sub = 0.9)
abline(v = log10(c(cz_small_thr, cz_large_thr)), lty = 2, col = "red")
text(x = log10(c(cz_small_thr, cz_large_thr)), y = max(max(hh$counts)), c(cz_small_thr, cz_large_thr), col = "darkred")
dev.off()

# remove large and small cells and create a new metacell object
metacell::mcell_mat_ignore_cells(
	new_mat_id = "nvec_filt",
	mat_id = "nvec",
	ig_cells=c(cells_small, cells_large))
# reload filtered matrix
mat = metacell::scdb_mat("nvec_filt")



#### Create metacell solution ####

# gene stats
# useful to select which genes to include in the clustering
metacell::mcell_add_gene_stat(
	gstat_id = "nvec_filt",
	mat_id = "nvec_filt",
	force = TRUE)


# define gene set to use for reclustering
# ignore genes that are not cell-type specific by definition, e.g. ribosomal genes? (might depend on your experiment)
genes_blacklisted = read.table("data/gene_list.ribosomal.txt")[,1]

# select gene set
metacell::mcell_gset_filter_multi(
	gstat_id = "nvec_filt",
	gset_id = "nvec_filt",
	T_tot = 100,
	T_top3 = 2,
	T_szcor = -0.05,
	T_niche = 0.01,
	force_new = TRUE,
	blacklist = genes_blacklisted)

# perform reclustering (to get a cgraph object)
metacell::mcell_add_cgraph_from_mat_bknn(
	mat_id = "nvec_filt",
	gset_id = "nvec_filt",
	graph_id = "nvec_filt",
	K = 100,
	dsamp = FALSE)

# coclustering graph via resampling (to get a coclust object)
metacell::mcell_coclust_from_graph_resamp(
	coc_id = "nvec_filt",
	graph_id = "nvec_filt",
	min_mc_size = 10,
	p_resamp = 0.75,
	n_resamp = 1000)

# balanced coclustering (to get a mc object)
metacell::mcell_mc_from_coclust_balanced(
	mc_id = "nvec_filt",
	coc_id = "nvec_filt",
	mat_id = "nvec_filt",
	K = 30,
	min_mc_size = 10,
	alpha = 2)

# load new metacells
mc = metacell::scdb_mc(id = "nvec_filt")

# are there any batch-correlated genes?
# find batch-correlated genes at the metacell level
batch_data_mc = sca_batch_correlated_genes_mc(mc_object = mc, mat_object = mat, cor_thr = 0.5, method = "spearman", do_plots = FALSE)
write.table(batch_data_mc$genes, sprintf("results_figs/%s.batch_genes.txt", "nvec"), row.names = FALSE, col.names = FALSE, quote = FALSE)



#### Order metacells ####

# use the metacell classification confusion matrix and dendrograms to order (and assign tentative colors) to metacells

# reload filtered objects
mc = metacell::scdb_mc(id = "nvec_filt")
mat = metacell::scdb_mat(id = "nvec_filt")

# obtain confusion matrix (mesuring cross-classification of metacells)
rec_confu_norm = mc_compute_norm_confu_matrix(
	mc_id = "nvec_filt", 
	graph_id = "nvec_filt", 
	max_deg = 100)

# dendrogram based on the confusion matrix, to decide how to group metacells
# visualise it before proceeding (see below)
rec_hc = mc_confusion_clustering(rec_confu_norm, clust_method = "ward.D2")
cut_rec_hc = min(as.dist(rec_hc)) + (max(as.dist(rec_hc)) - min(as.dist(rec_hc))) / 3 # an arbitrary heuristic to cut the tree, you can change it manually
mc_clusts = cutree(rec_hc, h=cut_rec_hc) # clusters after cutting the tree

# assign colors to each cluster of metacells
color_palette = colorRampPalette(c("magenta4","firebrick1","orange","khaki1","springgreen2","darkgreen","deepskyblue","cadetblue1","mediumblue","darkviolet","violet"))
cluster_colors = color_palette(n=length(table(mc_clusts[rec_hc$order])))
cluster_colors_per_mc = rep(cluster_colors, table(factor(mc_clusts[rec_hc$order],levels=unique(mc_clusts[rec_hc$order]))))
names(cluster_colors_per_mc) = 1:length(cluster_colors_per_mc)

# reorder metacells based on hierarchical clustering
mc_reord = metacell::mc_reorder(mc,rec_hc$order)
mc_reord@colors=cluster_colors_per_mc
# new metacell object, with proper order
metacell::scdb_add_mc(id = "nvec_filt_order", mc_reord)

# visualize dendrogram
# do you like the clusters?
pdf(sprintf("results_figs/%s.mc_dendrogram.pdf", "nvec"), height = 8, width = length(cluster_colors_per_mc) / 4)
rec_phy = ape::as.phylo(rec_hc)
ape::plot.phylo(rec_phy, direction="downwards", las=2, font=1, tip.color = "gray")
if ( length(unique(cutree(rec_hc, h = cut_rec_hc))) > 1  ) {
	rect.hclust(rec_hc,h=cut_rec_hc, border="darkgray")
}
text(x=1:length(cluster_colors_per_mc),y=0, names(cluster_colors_per_mc), col = cluster_colors_per_mc, srt = 270)
dev.off()

# plot confusion matrix (use the new reordered metacells)
# recalculate confusion matrix (with reordered metacells)
mc = metacell::scdb_mc(id = "nvec_filt_order")
rec_confu_norm = mc_compute_norm_confu_matrix(
	mc_id = "nvec_filt_order", 
	graph_id = "nvec_filt", 
	max_deg = 100)

pdf(sprintf("results_figs/%s.mc_confumatrix.pdf", "nvec"), width = 20, height = 20)
print(plot_complex_heatmap(
	sqrt(rec_confu_norm),
	name = "sqrt(norm confu)",
	color_min = 0, color_max = max(sqrt(rec_confu_norm)), 
	cluster_row = FALSE, cluster_col = FALSE,
	colors_row = mc@colors, colors_col = mc@colors,
	use_raster = FALSE))
dev.off()



#### Quality control ####

# plot useful stats about or new metacells: distribution of cell sizes, etc
mc = metacell::scdb_mc(id = "nvec_filt_order")
mat = metacell::scdb_mat(id = "nvec_filt")

# metacell size and counts
mc_stat = scp_plot_mc_size_counts(
	mc_object = mc, 
	mc_color = mc@colors, 
	mat_object = mat, 
	output_file = sprintf("results_figs/%s.mc_size_counts.pdf", "nvec"),
	return_table = TRUE)
write.table(mc_stat, sprintf("results_figs/%s.mc_size_counts.csv", "nvec"), quote = FALSE, sep = "\t", row.names = FALSE)



#### Plot metacell atlas ####

# create 2d projection and metacell-level and single cell-level heatmaps

# 2D PROJECTION
# create 2D projection of ordered metacells
metacell::mcell_mc2d_force_knn(
	mc2d_id = "nvec_filt",
	mc_id = "nvec_filt_order", 
	feats_gset = "nvec_filt", 
	graph_id = "nvec_filt")
# load mc2d object
mc2d = metacell::scdb_mc2d("nvec_filt")

# plot 2D projection of ordered metacells
scp_plot_sc_2d(
	mc2d = mc2d,
	mc = mc,
	plot_mcs=TRUE,
	plot_mc_name=TRUE,
	plot_edges=TRUE,
	output_file=sprintf("results_figs/%s.2dproj.pdf",  "nvec"),
	width=10, height=10)


# gene expression heatmaps
# first, load some annotations for your genes (rownames are genes)
gene_annotations = read.table("data/gene_names.tsv", sep = "\t", col.names = c("gene","gene_name","domains"), row.names = 1)

# select some marker genes that are metacell-specific (based on FC>1.5), up to 20 genes per cluster
marker_data_list = scp_plot_cmod_markers_select(
	mc_object = mc,
	gene_annot_file = gene_annotations,
	per_clust_genes = 20,
	gene_min_fold = 2,
	gene_font_size = 5,
	clust_ord = 1:ncol(mc@mc_fp))

# which genes are TFs? They'll be highlighited in the heatmap
gene_annotations_tfs = read.table("data/gene_list.TFs.tsv", row.names = 1)
gene_list_tfs = rownames(gene_annotations_tfs)

# draw metacell heatmap
scp_plot_cmod_markers_mc(
	marker_data_list = marker_data_list,
	output_file = sprintf("results_figs/%s.exp_global_mc.pdf", "nvec"),
	heatmap_colors = c("white","orange","orangered2","#520c52"),
	clust_col = mc@colors,
	width = max(10, round(length(marker_data_list$clust_ord) / 10)), # these width and height heuristics are designed to look good in most cases
	height = round(length(marker_data_list$genes) / 15) + 1,
	show_gene_names = TRUE,
	highlight_genes = gene_list_tfs,
	gene_chr_limit = 60,
	show_clust_borders = TRUE,
	use_raster=TRUE,
	max_expression_fc = 3,
	min_expression_fc = 1
)

# a similar heatmap, at the single-cell level
scp_plot_cmod_markers_sc(
	marker_data_list = marker_data_list,
	mc_object = mc,
	mat_object = mat,
	output_file = sprintf("results_figs//%s.exp_global_sc.pdf",  "nvec"),
	heatmap_colors = c("white","orange","orangered2","#520c52"),
	clust_col = mc@colors,
	width = max(12, round(length(marker_data_list$clust_ord) / 2)), # these width and height heuristics are designed to look good in most cases
	height = round(length(marker_data_list$genes) / 15) + 1,
	show_gene_names = TRUE,
	highlight_genes = gene_list_tfs,
	smoothen = 5,
	gene_chr_limit = 60,
	show_clust_borders=TRUE,
	use_raster=TRUE,
	max_expression_fc = 3,
	min_expression_fc = 1
)

# we could also plot the expression of all tfs
mc_counts = sca_mc_gene_counts(mc, mat, 0)  # counts of each gene per metacell
gene_subset_m = scp_barplot_heatmap_markers(
	mc_object = mc,
	mat_object = mat,
	mc_counts = mc_counts,
	markers_file = gene_annotations_tfs,
	heatmap_colors =  c("white","orange","orangered2","#520c52"),
	output_file_heatmap = sprintf("results_figs/%s.markers_TFs.heatmap.pdf", "nvec"),
	output_file_barplot = sprintf("results_figs/%s.markers_TFs.barplot.pdf", "nvec"),
	T_totumi = 10,
	width = 10,
	height = NULL,
	use_raster = TRUE,
	min_gene_fc = 1.5,
	min_expression_fc = 1,
	max_expression_fc = 3,
	mc_color = mc@colors,
	print_barplots = TRUE
)

# save our tentative metacell annotation table
mc_t = data.frame(
	metacell = colnames(mc@mc_fp),
	cell_type = paste("cell_type_", as.numeric(factor(mc@colors, levels = unique(mc@colors))) , sep = ""), # this will name tentative cell types according to their assigned colors
	color = mc@colors
)
write.table(mc_t, file = "data/annotation_mc.it0.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# you can now start annotating your metacells and metacell clusters!