library(Seurat)
library(harmony)
library(magrittr)
library(ggplot2)
library(stringr)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12")
dir.create(TARGET_dir)

object = readRDS(file.path("results", ANALYSIS_VERSION, "early_wt12_integrated/BDGP_automated_annotation_object.rds"))

manual_tab = read.csv(file.path(TARGET_dir, 'manualCellType_3.csv'))
object@meta.data$manual_celltypes = NULL

# dorsal epidermis
VlnPlot(object, features = 'yellow-e3', pt.size = 0)
FeaturePlot(object, features = 'yellow-e3')

VlnPlot(object, features = 'ft', pt.size = 0)
FeaturePlot(object, features = 'ft')

VlnPlot(object, features = 'CG7296', pt.size = 0)
FeaturePlot(object, features = 'CG7296')

for(temp_cluster in unique(manual_tab$cluster)) { 
  object@meta.data[object@meta.data$seurat_clusters == temp_cluster, 'manual_celltypes'] = trimws(manual_tab[manual_tab$cluster == temp_cluster, 'annotation'])
}

p = DimPlot(object, group.by = 'manual_celltypes', label = TRUE)
ggsave(filename = file.path(TARGET_dir, "manual_celltypes.png"), plot = p, width = 10, height = 8)
saveRDS(object, file = file.path(TARGET_dir, "manual_celltype_object1.rds"))

