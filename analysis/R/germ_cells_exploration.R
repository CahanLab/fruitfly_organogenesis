library(Seurat)
library(enrichR)
enrichR::setEnrichrSite("FlyEnrichr")

TARGET_dir = file.path("results", ANALYSIS_VERSION, "germ_cell_exploration")
dir.create(TARGET_dir)

early_wt_object = readRDS("results/v18/manual_annotation_early_wt12/manual_celltype_object1.rds")
sub_early_object = early_wt_object[, early_wt_object$seurat_clusters %in% c(32, 24)]

DimPlot(sub_early_object)
markers = SeuratWrappers::RunPresto(sub_early_object, ident.1 = 32, logfc.threshold = 0, min.pct = 0.1, group.by = 'seurat_clusters') 

markers_32 = markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0, ]
markers_24 = markers[markers$p_val_adj < 0.05 & markers$avg_log2FC < 0, ]

write.csv(markers_32, file = file.path(TARGET_dir, "markers_32_genes.csv"))
write.csv(markers_24, file = file.path(TARGET_dir, "markers_24_genes.csv"))

enrichment_results = enrichR::enrichr(
  genes = rownames(markers_24), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018"
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
write.csv(biological_analysis, file = file.path(TARGET_dir, "GO_24_enrichment.csv"))

# cluster 32 
enrichment_results = enrichR::enrichr(
  genes = rownames(markers_32), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018"
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
write.csv(biological_analysis, file = file.path(TARGET_dir, "GO_32_enrichment.csv"))

#######################################
late_wt_object = readRDS("results/v18/manual_annotation_wt13/manual_celltype_object3.rds")
cds = readRDS("results/v18/refined_wt_late_early_germ/monocle3_no_batch_correct_object.rds")
late_wt_object@meta.data$weird_cluster = 'No'
late_wt_object@meta.data[names(monocle3::clusters(cds))[monocle3::clusters(cds) == 1], 'weird_cluster'] = 'Yes'
DimPlot(late_wt_object, group.by = 'weird_cluster')
ggsave(filename = file.path(TARGET_dir, 'weird_germ_cluster.png'), width = 8, height = 6)

sub_late_object = subset(late_wt_object, subset = manual_celltypes == 'Germ Cell')
DimPlot(sub_late_object, group.by = 'manual_celltypes')
old_sub_late_object = sub_late_object

sub_late_object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
sub_late_object %<>% ScaleData(features = VariableFeatures(object = sub_late_object))
sub_late_object %<>% RunPCA(features = VariableFeatures(object = sub_late_object), npcs = 100)
num_pc = 20
sub_late_object %<>% FindNeighbors(dims = 1:num_pc)
sub_late_object %<>% RunUMAP(dim = 1:num_pc)
sub_late_object %<>% FindClusters(resolution = 0.1)
DimPlot(sub_late_object)
old_sub_late_object@meta.data$new_cluster = sub_late_object@meta.data$seurat_clusters

DimPlot(old_sub_late_object, group.by = 'new_cluster')
ggsave(filename = file.path(TARGET_dir, 'new_old_late_cluster.png'), width = 8, height = 6)

markers = SeuratWrappers::RunPresto(old_sub_late_object, ident.1 = 1, logfc.threshold = 0, min.pct = 0.1, group.by = 'new_cluster') 

markers_1 = markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0, ]
markers_2 = markers[markers$p_val_adj < 0.05 & markers$avg_log2FC < 0, ]

write.csv(markers_1, file = file.path(TARGET_dir, "late_markers_1_genes.csv"))
write.csv(markers_2, file = file.path(TARGET_dir, "late_markers_2_genes.csv"))

enrichment_results = enrichR::enrichr(
  genes = rownames(markers_1), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018"
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
write.csv(biological_analysis, file = file.path(TARGET_dir, "GO_late_1_enrichment.csv"))

# cluster 32 
enrichment_results = enrichR::enrichr(
  genes = rownames(markers_2), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018"
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
write.csv(biological_analysis, file = file.path(TARGET_dir, "GO_late_2_enrichment.csv"))





