
TARGET_dir = file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13")

##### first round of cell typing #####
object = readRDS(file.path("results", ANALYSIS_VERSION, "wt13_integrated/BDGP_automated_annotation_object.rds"))

manual_tab = read.csv(file.path('accessory_data/manual_annotations_v18/wt13/', 'manualCellType.csv'))

object@meta.data$manual_celltypes = NULL

for(temp_cluster in unique(manual_tab$cluster)) { 
  object@meta.data[object@meta.data$seurat_clusters == temp_cluster, 'manual_celltypes'] = trimws(manual_tab[manual_tab$cluster == temp_cluster, 'annotation'])
}

marker_genes = read.csv(file.path(TARGET_dir, 'manualMarkerGenes.csv'))
marker_genes_list = list()

for(unique_ct in unique(marker_genes$annotation)) {
  sub_marker_genes = marker_genes[marker_genes$annotation == unique_ct, ]
  marker_genes_list[[trimws(unique_ct)]] = intersect(as.vector(sub_marker_genes$marker_genes), rownames(object@assays$RNA@data))
  object = AddModuleScore(
    object = object,
    features = list(marker_genes_list[[trimws(unique_ct)]]),
    name = paste0(stringr::str_replace_all(trimws(unique_ct), pattern = " ", replacement = "_"), '_ModuleScore')
  )
}

FeaturePlot(object, features = 'Glia_ModuleScore1')
VlnPlot(object, features = 'Glia_Midline_ModuleScore1', pt.size = 0)
object@meta.data[object@meta.data$seurat_clusters == 26, 'manual_celltypes'] = "Muscle System"
p = DimPlot(object, group.by = 'manual_celltypes', label = TRUE)
ggsave(filename = file.path(TARGET_dir, 'cell_type_UMAP.png'), plot = p, width = 18, height = 9)
saveRDS(object, file = file.path(TARGET_dir, 'manual_celltype_object.rds'))

##### second round of cell typing #####
object = readRDS(file.path(TARGET_dir, 'manual_celltype_object.rds'))
Seurat::DimPlot(object, group.by = 'manual_celltypes')

# rename all the epidermis subtype as epidermis 
object@meta.data[object@meta.data$manual_celltypes == 'Esophagus', 'manual_celltypes'] = 'Epidermis'
object@meta.data[object@meta.data$manual_celltypes == 'Ventral Epidermis', 'manual_celltypes'] = 'Epidermis'
object@meta.data[object@meta.data$manual_celltypes == 'Hypopharynx', 'manual_celltypes'] = 'Epidermis'
object@meta.data[object@meta.data$manual_celltypes == 'Dorsal Epidermis', 'manual_celltypes'] = 'Epidermis'
object@meta.data[object@meta.data$manual_celltypes == 'Epipharynx', 'manual_celltypes'] = 'Epidermis'

p = DimPlot(object, group.by = 'manual_celltypes', label = TRUE)
ggsave(filename = file.path(TARGET_dir, 'cell_type_UMAP2.png'), plot = p, width = 18, height = 9)

saveRDS(object, file.path(TARGET_dir, 'manual_celltype_object2.rds'))

##### third round of cell typing ##### 
object = readRDS(file.path(TARGET_dir, "manual_celltype_object2.rds"))
object@meta.data[object@meta.data$seurat_clusters == 26, 'manual_celltypes'] = 'Hindgut Muscle'
object@meta.data[object@meta.data$seurat_clusters == 13, 'manual_celltypes'] = 'Pharyngeal Muscle'
p = DimPlot(object, group.by = 'manual_celltypes', label = TRUE)
ggsave(filename = file.path(TARGET_dir, 'cell_type_UMAP3.png'), plot = p, width = 18, height = 9)
saveRDS(object, file.path(TARGET_dir, 'manual_celltype_object3.rds'))

##### fourth round of cell typing #####
object = readRDS(file.path(TARGET_dir, 'manual_celltype_object.rds'))
object@meta.data[object@meta.data$seurat_clusters == 26, 'manual_celltypes'] = 'Hindgut Muscle'
object@meta.data[object@meta.data$seurat_clusters == 13, 'manual_celltypes'] = 'Pharyngeal Muscle'
object@meta.data[object@meta.data$manual_celltypes == 'Germ Cell', 'manual_celltypes'] = 'Germ Cells'
object@meta.data[object@meta.data$manual_celltypes == 'Longitudinal/Caudal Visceral Mesoderm', 'manual_celltypes'] = 'Longitudinal Visceral Mesoderm'
object@meta.data[object@meta.data$manual_celltypes == 'Visceral Mesoderm', 'manual_celltypes'] = 'Circular Visceral Mesoderm'
object@meta.data[object@meta.data$seurat_clusters == 15 | object@meta.data$seurat_clusters == 4, 'manual_celltypes'] = 'Early CNS'
object@meta.data[object@meta.data$seurat_clusters %in% c(33, 3, 11), 'manual_celltypes'] = 'Late CNS'
DimPlot(object, group.by = 'manual_celltypes', label = TRUE)
saveRDS(object, file.path(TARGET_dir, 'manual_celltype_object4.rds'))


