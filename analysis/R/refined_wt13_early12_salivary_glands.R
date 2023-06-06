# analyze the salivary gland cells between early and late 
set.seed(123)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland")
dir.create(TARGET_dir)

##### load in the latest version of the files #####
wt_late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object4.rds"))
sub_wt_late_object = subset(wt_late_object, subset = manual_celltypes == "Salivary Gland")
sub_wt_late_object$experimental_condition = 'late'

wt_early_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))
sub_wt_early_object = subset(wt_early_object, subset = manual_celltypes == 'Salivary Gland')
sub_wt_early_object$experimental_condition = 'early'

##### compile the early and late data using seurat #####
combined_ct_object = merge(sub_wt_early_object, sub_wt_late_object)
object = Seurat::CreateSeuratObject(combined_ct_object@assays$RNA@counts, project = 'salivary_glands')
object@meta.data$experimental_condition = combined_ct_object@meta.data$experimental_condition
object@meta.data$batch = paste0(combined_ct_object@meta.data$experimental_condition, "_", combined_ct_object@meta.data$batch)
object = Seurat::NormalizeData(object)
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

##### convert the seurat data into monocle3 data and preprocess #####
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

# plot out an early SG gene 
plot_cells(cds,
           genes=c("eyg"),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_eyg.png'), width = 8, height = 6)

# cluster the cells and performed DE genes analysis using presto 
cds <- cluster_cells(cds, resolution = 1e-2)
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1)

rank_sum_results = presto::wilcoxauc(normalized_counts(cds), cds@clusters$UMAP$clusters)
write.csv(rank_sum_results, file = file.path(TARGET_dir, "rank_sum_test.csv"))

# plot out the phase 
plot_cells(cds, color_cells_by = "Phase", label_cell_groups = FALSE, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_phase.png'), width = 8, height = 6)

# pseudotime ordering using monocle3 
cds <- learn_graph(cds, use_partition = FALSE)
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
saveRDS(cds, file = file.path(TARGET_dir, "monocle3_no_batch_correct_object.rds"))

# plot out some genes for sanity checks 
plot_cells(cds,
           genes=c("eyg", "CrebA", 'rib', 'wbl', 'toe', 'sage', 'fkh', 'sano', 'pip', 'RpS17', 'RpL24', 'RpL41', 'RpL13A'),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)

##### perform GSEA on the two clusters #####
pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')
rank_sum_test = read.csv(file.path(TARGET_dir, 'rank_sum_test.csv'), row.names = 1)
rank_sum_test = rank_sum_test[rank_sum_test$pct_in > 10 | rank_sum_test$pct_out > 10, ]
sub_rank_sum_test = rank_sum_test[rank_sum_test$group == 1, ]
ranks <- sub_rank_sum_test$logFC
names(ranks) <- sub_rank_sum_test$feature
fgseaRes <- fgsea(pathways = pathway_list, 
                  stats = ranks,
                  minSize=10,
                  maxSize=500)

fgseaRes = data.frame(fgseaRes) 
fgseaRes = apply(fgseaRes,2,as.character)
fgseaRes = as.data.frame(fgseaRes)
fgseaRes = fgseaRes[!is.na(fgseaRes$padj), ]
fgseaRes$NES = as.numeric(fgseaRes$NES)
fgseaRes = fgseaRes[fgseaRes$NES > 0, ]

write.csv(fgseaRes, file = file.path(TARGET_dir, 'late_gsea_results.csv'))

sub_rank_sum_test = rank_sum_test[rank_sum_test$group == 2, ]
ranks <- sub_rank_sum_test$logFC
names(ranks) <- sub_rank_sum_test$feature
fgseaRes <- fgsea(pathways = pathway_list, 
                  stats = ranks,
                  minSize=10,
                  maxSize=500)

fgseaRes = data.frame(fgseaRes)
fgseaRes = apply(fgseaRes,2,as.character)
fgseaRes = as.data.frame(fgseaRes)
fgseaRes = fgseaRes[!is.na(fgseaRes$padj), ]
fgseaRes$NES = as.numeric(fgseaRes$NES)
fgseaRes = fgseaRes[fgseaRes$NES > 0, ]

write.csv(fgseaRes, file = file.path(TARGET_dir, 'early_gsea_results.csv'))
