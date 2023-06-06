# analysis of trachea cells 

set.seed(123)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea")
dir.create(TARGET_dir)

##### load in the appropriate data #####
wt_late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object4.rds"))
sub_wt_late_object = subset(wt_late_object, subset = manual_celltypes == "Trachea")
sub_wt_late_object$experimental_condition = 'late'

wt_early_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))
sub_wt_early_object = subset(wt_early_object, subset = manual_celltypes == 'Trachea')
sub_wt_early_object$experimental_condition = 'early'

combined_ct_object = merge(sub_wt_early_object, sub_wt_late_object)
object = Seurat::CreateSeuratObject(combined_ct_object@assays$RNA@counts, project = 'Trachea')
object@meta.data$experimental_condition = combined_ct_object@meta.data$experimental_condition
object@meta.data$batch = paste0(combined_ct_object@meta.data$experimental_condition, "_", combined_ct_object@meta.data$batch)

object = Seurat::NormalizeData(object)
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

##### convert to monocle3 and preprocess #####
expression_matrix = object@assays$RNA@counts
cell_metadata = object@meta.data
gene_annotation = object@assays$RNA@meta.features
gene_annotation$gene_short_name = rownames(gene_annotation)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "batch", cell_size = 1, label_cell_groups = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_UMAP.png'), width = 8, height = 6)

# cluster the cells 
cds <- cluster_cells(cds, resolution = 5e-3)
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_cluster.png'), width = 8, height = 6)

# calculate DE genes using preseto package
rank_sum_results = presto::wilcoxauc(normalized_counts(cds), cds@clusters$UMAP$clusters)
write.csv(rank_sum_results, file = file.path(TARGET_dir, "rank_sum_test.csv"))

plot_cells(cds, color_cells_by = "Phase", label_cell_groups = FALSE, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_phase.png'), width = 8, height = 6)

cds <- learn_graph(cds, use_partition = FALSE, learn_graph_control = list("minimal_branch_len" = 18))

plot_cells(cds,
           color_cells_by = "batch",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_pt.png'), width = 8, height = 6)

plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1, cell_size = 1, show_trajectory_graph = FALSE)


cds@colData$cell_type = NA
cds@colData$cluster_id = monocle3::clusters(cds)
cds@colData[cds@colData$cluster_id %in% c(3, 7, 6), 'cell_type'] = 'Early Trachea Cells'
cds@colData[cds@colData$cluster_id %in% c(5, 4, 2), 'cell_type'] = 'Middle Trachea Cells'
cds@colData[cds@colData$cluster_id %in% c(8, 1), 'cell_type'] = 'Late Trachea Cells'
cds@colData[cds@colData$cluster_id %in% c(9), 'cell_type'] = 'Branching Trachea Cells'

saveRDS(cds, file = file.path(TARGET_dir, "monocle3_no_batch_correct_object.rds"))

##### perform GSEA on the sub-populations #####
pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')
rank_sum_test = read.csv(file.path(TARGET_dir, 'rank_sum_test.csv'), row.names = 1)
rank_sum_test = rank_sum_test[rank_sum_test$pct_in > 10 | rank_sum_test$pct_out > 10, ]
for(ct in unique(rank_sum_test$group)) {
  sub_rank_sum_test = rank_sum_test[rank_sum_test$group == ct, ]
  ranks <- sub_rank_sum_test$logFC
  names(ranks) <- sub_rank_sum_test$feature
  fgseaRes <- fgsea(pathways = pathway_list, 
                    stats = ranks,
                    minSize=10,
                    maxSize=500)
  
  fgseaRes = data.frame(fgseaRes)
  fgseaRes = apply(fgseaRes,2,as.character)
  fgseaRes = as.data.frame(fgseaRes)
  fgseaRes$padj = as.numeric(fgseaRes$padj)
  fgseaRes = fgseaRes[!is.na(fgseaRes$padj), ]
  fgseaRes$NES = as.numeric(fgseaRes$NES)
  fgseaRes = fgseaRes[fgseaRes$NES > 0, ]
  write.csv(fgseaRes, file = file.path(TARGET_dir, paste0(ct, '_gsea_results.csv')))
}

