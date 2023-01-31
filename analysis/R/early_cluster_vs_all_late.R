library(Seurat)
library(SeuratWrappers)
library(presto)
library(fgsea)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "early_celltype_vs_late_all")
dir.create(TARGET_dir)

early_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12", "manual_celltype_object1.rds"))
late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13", "manual_celltype_object3.rds"))
late_object@meta.data$target_ct = 'late'

early_object@meta.data$fine_celltype = paste0(early_object@meta.data$seurat_clusters, "_", early_object@meta.data$manual_celltypes)

pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')

for(fine_cluster in unique(early_object@meta.data$fine_celltype)) {
  print(fine_cluster)
  sub_early_object = subset(early_object, subset = fine_celltype == fine_cluster)
  sub_early_object@meta.data$target_ct = fine_cluster
  combined_object = merge(sub_early_object, late_object)
  object = Seurat::CreateSeuratObject(combined_object@assays$RNA@counts, project = 'early_ct_vs_late_all')
  object@meta.data$experimental_condition = combined_object@meta.data$target_ct
  object = Seurat::NormalizeData(object)
  markers = SeuratWrappers::RunPresto(object, ident.1 = fine_cluster, logfc.threshold = 0, min.pct = 0.1, group.by = 'experimental_condition') 
  write.csv(x = markers, file = file.path(TARGET_dir, paste0(fine_cluster, "_early_vs_late_all.csv")))
}

dir.create(file.path(TARGET_dir, "GSEA_results"))
for(fine_cluster in unique(early_object@meta.data$fine_celltype)) {
  print(fine_cluster)
  markers = read.csv(file.path(TARGET_dir, paste0(fine_cluster, "_early_vs_late_all.csv")), row.names = 1)
  
  ranks <- markers$avg_log2FC
  names(ranks) <- rownames(markers)
  
  fgseaRes <- fgsea(pathways = pathway_list, 
                    stats = ranks,
                    minSize=10,
                    maxSize=500)
  fgseaRes = data.frame(fgseaRes)
  fgseaRes <- apply(fgseaRes,2,as.character)
  
  withr::with_dir(file.path(TARGET_dir, "GSEA_results"), {
    write.csv(fgseaRes, file = paste0(fine_cluster, '_gsea_results_wt.csv'))
  })
}
