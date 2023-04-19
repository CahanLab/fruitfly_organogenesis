library(Seurat)
library(viridis)

Seroka_object = readRDS("accessory_data/Doe_Drosophila_Embryo_Atlas/script/curated_embryo_Doe.rds")
DimPlot(Seroka_object)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "seroka_sg_genes")
dir.create(TARGET_dir)

FeaturePlot(Seroka_object, features = 'pip')
ggsave(file.path(TARGET_dir, 'pip_expression.png'))
VlnPlot(Seroka_object, features = 'pip')

Seroka_object@meta.data$sg_ct = 'other'
Seroka_object@meta.data[Seroka_object@meta.data$seurat_clusters == '52', 'sg_ct'] = 'sg_1'
Seroka_object@meta.data[Seroka_object@meta.data$seurat_clusters == '60', 'sg_ct'] = 'sg_2'

DimPlot(Seroka_object, group.by = 'sg_ct')
ggsave(file.path(TARGET_dir, 'new_labels.png'))

VlnPlot(Seroka_object, features = 'pip', group.by = 'sg_ct')

markers = SeuratWrappers::RunPresto(Seroka_object, ident.1 = 'sg_1', logfc.threshold = 0, min.pct = 0.1, group.by = 'sg_ct') 
write.csv(markers, file = file.path(TARGET_dir, 'sg_1_marker.csv'))

markers = SeuratWrappers::RunPresto(Seroka_object, ident.1 = 'sg_2', logfc.threshold = 0, min.pct = 0.1, group.by = 'sg_ct') 
write.csv(markers, file = file.path(TARGET_dir, 'sg_2_marker.csv'))
