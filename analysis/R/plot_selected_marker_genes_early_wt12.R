# make plots of rare cell types in stage 10-12 
# Supp Fig 13-16

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'wt_early_selected_marker_genes')
dir.create(TARGET_dir, recursive = TRUE)

object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))

marker_gene_list = list()
marker_gene_list[['Crystal Cells']] = c('PPO1', 'PPO2', 'CG17109')
marker_gene_list[['Malpighian tubules']] = c('CG3699', 'Fst')
marker_gene_list[['Apodemes']] = c('CG7296', 'TwdlM', 'TwdlB', 'CG7299', 'TwdlL')
marker_gene_list[['Hindgut']] = c('otp', 'CG17780', 'CG13082')
marker_gene_list[['Salivary Gland']] = c('sano', 'pip', 'CG9005', 'PH4alphaSG2')

withr::with_dir(TARGET_dir, {
  for(ct in names(marker_gene_list)) {
    for(gene in marker_gene_list[[ct]]) {
      if(gene %in% rownames(object)) {
        p = Seurat::VlnPlot(object, features = gene, group.by = 'manual_celltypes', pt.size = 0) + 
          ylab("Normalized Expression") + 
          xlab("Stage 10 - 12: Drosophila Embryonic Cell Types") + 
          ggtitle(paste0(gene, " (", ct, ")")) +
          theme(text = element_text(size = 16), legend.text=element_text(size=24), legend.position = "none")
        ggsave(filename = paste0("vlnplot_", ct, "_", gene, ".png"), plot = p, width = 17, height = 10)
        
        p = Seurat::FeaturePlot(object, features = gene) + 
          ggtitle(paste0(gene, " (", ct, ")")) + 
          theme(text = element_text(size = 16), legend.text=element_text(size=24), legend.position = 'none')
        ggsave(filename = paste0("UMAP_", ct, "_", gene, ".png"), plot = p, width = 10, height = 10)
        
      }
    }
  }
})

