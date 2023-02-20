library(ggplot2)
library(RColorBrewer)
library(ggdendroplot)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'wt_late_selected_marker_genes')
dir.create(TARGET_dir, recursive = TRUE)

object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object4.rds"))

marker_gene_list = list()
marker_gene_list[['Dorsal Vessel']] = c('tin', 'prc', 'tup', 'Him')
marker_gene_list[['Garland Cells']] = c('CG42255', 'beta4GalNacTA', 'Amnionless', 'CG15209')
marker_gene_list[['Esophagus']] = c('ect', 'amd', 'CG13631', 'Osi7', 'Osi14')
marker_gene_list[['Hypopharynx']] = c('CG6118', 'bbg', 'CG10211')
marker_gene_list[['Epipharynx']] = c('yellow-e3', 'optix', 'sprt')
marker_gene_list[['Hindgut Muscle']] = c('wnt4', 'AbdB', 'bin', 'up')
marker_gene_list[['Pharyngeal Muscle']] = c('bt', 'up', 'CG5080')
marker_gene_list[['Optic Lobe']] = c('E(Spl)M5-HLH', 'Obp99a', 'SoxN')
marker_gene_list[['Malpighian tubules']] = c('CG31272', 'Bw', 'ZnT35C')


withr::with_dir(TARGET_dir, {
  for(ct in names(marker_gene_list)) {
    for(gene in marker_gene_list[[ct]]) {
      if(gene %in% rownames(object)) {
        p = Seurat::VlnPlot(object, features = gene, group.by = 'manual_celltypes', pt.size = 0) + 
          ylab("Normalized Expression") + 
          xlab("Stage 13 - 16: Drosophila Embryonic Cell Types") + 
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

