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

##### plot out the violin plot of pseudotime without the extra unknown cells #####
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ", "monocle3_no_batch_correct_object.rds"))
UMAP_coord = cds@int_colData$reducedDims$UMAP

colnames(UMAP_coord) = c("UMAP_1", "UMAP_2")
UMAP_coord = as.data.frame(UMAP_coord)
UMAP_coord$clusters = as.vector(monocle3::clusters(cds))
UMAP_coord$pseudotime = as.vector(monocle3::pseudotime(cds))
UMAP_coord$batch = as.vector(cds@colData$batch)

UMAP_coord[UMAP_coord$batch == 'early_rep_1', 'batch'] = 'Stage 10-12 rep 1'
UMAP_coord[UMAP_coord$batch == 'early_rep_2', 'batch'] = 'Stage 10-12 rep 2'
UMAP_coord[UMAP_coord$batch == 'late_rep_1', 'batch'] = 'Stage 13-16 rep 1'
UMAP_coord[UMAP_coord$batch == 'late_rep_3', 'batch'] = 'Stage 13-16 rep 2'

UMAP_coord$cell_type = UMAP_coord$clusters
UMAP_coord[UMAP_coord$clusters == '1', 'cell_type'] = 'Unknown 1'
UMAP_coord[UMAP_coord$clusters == '3', 'cell_type'] = 'Unknown 2'
UMAP_coord[UMAP_coord$clusters == '5', 'cell_type'] = 'Early Germ Cells'
UMAP_coord[UMAP_coord$clusters == '6', 'cell_type'] = 'Interm. Germ Cells 1'
UMAP_coord[UMAP_coord$clusters == '2', 'cell_type'] = 'Interm. Germ Cells 2'
UMAP_coord[UMAP_coord$clusters == '4', 'cell_type'] = 'Late Germ Cells'

UMAP_coord = UMAP_coord[UMAP_coord$cell_type != 'Unknown 1', ]
UMAP_coord = UMAP_coord[UMAP_coord$cell_type != 'Unknown 2', ]

p = ggplot(UMAP_coord, aes(x=reorder(batch, pseudotime), y=pseudotime, fill = batch)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set1') + 
  ylab("Pseudotime") + 
  xlab("Batch") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_pseudotime_germ_no_unknown.png"), plot = p, width = 8, height = 6)

##### look at the TFs in the SG on UMAP #####
cds = readRDS("results/v18/figure_plots/clean_sharable_data/SG_specific/SG_monocle3_object.rds")
interesting_genes = c('toe', 'hkb', 'eyg', 'Scr', 'trh', 'Dr', 'sens', 'brk', 'D', 'bowl','fkh','sage', 'SoxN', 'ttk', 'CrebA','rib')
p = plot_cells(cds, genes = interesting_genes, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(filename = file.path(TARGET_dir, "TFs_UMAP_SG.png"))
