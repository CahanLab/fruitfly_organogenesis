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

##### verify if there is a difference between early and late tracheal tip cells #####
library(cowplot)

cds = readRDS("results/v18/figure_plots/clean_sharable_data/trachea_specific/trachea_monocle3_object.rds")
norm_exp = monocle3::normalized_counts(cds)
meta_tab = cds@colData
meta_tab = meta_tab[meta_tab$subtypes == 'Tracheal Tip Cells', ]
norm_exp = norm_exp[, rownames(meta_tab)]
rank_sum_results = presto::wilcoxauc(norm_exp, meta_tab$experimental_condition)
write.csv(rank_sum_results, file = file.path(TARGET_dir, 'rank_sum_results_early_late_tip.csv'))
early_sum_results = rank_sum_results[rank_sum_results$group == 'early', ]
early_sum_results$abs_logFC = abs(early_sum_results$logFC)
early_sum_results$category = NA
early_sum_results[early_sum_results$logFC > 0, 'category'] = 'stage10-12 Tip Cells'
early_sum_results[early_sum_results$logFC < 0, 'category'] = 'stage13-16 Tip Cells'
early_sum_results = early_sum_results[order(early_sum_results$abs_logFC, decreasing = TRUE), ]
early_sum_results = early_sum_results[1:20, ]
p = ggplot(data = early_sum_results, aes(x = reorder(feature, abs_logFC), y = abs_logFC, fill = category)) + 
  geom_bar(stat = 'identity') + coord_flip() + 
  ylab("absolute logFC") + 
  xlab("differentially expressed genes") + 
  theme_half_open()
ggsave(filename = file.path(TARGET_dir, "DE_genes_tip_cells.png"), plot = p, width = 5, height = 5)

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



