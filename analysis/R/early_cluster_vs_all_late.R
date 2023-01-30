library(Seurat)
library(SeuratWrappers)
library(presto)

early_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12", "manual_celltype_object1.rds"))
late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13", "manual_celltype_object3.rds"))

early_object@meta.data$fine_celltype = paste0(early_object@meta.data$seurat_clusters, "_", early_object@meta.data$manual_celltypes)

for(fine_cluster in unique(early_object@meta.data$fine_celltype)) {
  print(fine_cluster)
}
