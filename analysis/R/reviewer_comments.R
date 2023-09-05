# this script contain analysis to address reviewer comments 
# but the end results may not end up in the manuscript 
set.seed(123)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "reviewer_comments")
dir.create(TARGET_dir)

##### look at SG UMAP in 3D #####
cds = readRDS("results/v18/figure_plots/clean_sharable_data/SG_specific/SG_monocle3_object.rds")
cds@colData[cds@colData$batch == 'late_rep_3', 'batch'] = 'late_rep_2'
cds <- monocle3::reduce_dimension(cds, max_components = 3)
plot_cells_3d(cds, color_cells_by = "batch")

plot_cells(cds, x = 1, y = 3, cell_size = 0.9, group_label_size = 4, color_cells_by = 'batch')
ggsave(file.path(TARGET_dir, "3D_UMAP1.png"), width = 5, height = 5)

plot_cells(cds, x = 2, y = 3, cell_size = 0.9, group_label_size = 4, color_cells_by = 'batch')
ggsave(file.path(TARGET_dir, "3D_UMAP2.png"), width = 5, height = 5)
