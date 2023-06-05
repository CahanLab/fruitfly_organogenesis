# plot out nfeatures and log10_ncounts across batches 
TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'wt_late_statistics')
dir.create(TARGET_dir, recursive = TRUE)

object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object4.rds"))
object@meta.data[object@meta.data$batch == 'rep_3', 'batch'] = 'rep_2'

withr::with_dir(TARGET_dir, {
  p = Seurat::VlnPlot(object, features = 'nFeature_RNA', group.by = 'batch', pt.size = 0) + 
      ylab("Number of Genes Expressed") + 
      xlab("Batch") + 
    theme(text = element_text(size = 22), legend.text=element_text(size=22))
  ggsave(filename = 'nFeatures_RNA.png', plot = p, width = 8, height = 6)
  
  p = Seurat::VlnPlot(object, features = 'log10_nCount_RNA', group.by = 'batch', pt.size = 0) + 
      ylab("log10(# of Counts)") + 
      xlab("Batch") + 
    theme(text = element_text(size = 22), legend.text=element_text(size=22))
  ggsave(filename = 'log10_nCount_RNA.png', plot = p, width = 8, height = 6)
})

#> table(object@meta.data$batch)
#rep_1 rep_2 
#22054 20663 


