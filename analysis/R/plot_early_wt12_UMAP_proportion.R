# plot out UMAP and cell type proportion 
# Fig 1D, 1E
# Supp Fig 1C, 1D
TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'wt_early_UMAP_proportion')
dir.create(TARGET_dir, recursive = TRUE)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
qual_col_pals = qual_col_pals[rownames(qual_col_pals) != 'Dark2', ]
qual_col_pals = qual_col_pals[rownames(qual_col_pals) != 'Set2', ]

col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = col_vector[col_vector != col_vector[4]]
col_vector = col_vector[col_vector != col_vector[4]]

wt_object = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_early_wt12/manual_celltype_object1.rds'))

# plot the UMAP with cell type labels, seurat clusters and batch information
withr::with_dir(TARGET_dir, {
  p = DimPlot(wt_object, group.by = 'manual_celltypes', label = FALSE, label.size = 5) +
    ggtitle("Stage 10 – 12 Drosophila Embryonic Cell Type Labels") + 
    theme(text = element_text(size = 18))
  ggsave(file.path("cell_type_UMAP_Unlabelled.png"), plot = p, width = 17, height = 10)
  
  p = DimPlot(wt_object, group.by = 'manual_celltypes', label = TRUE, label.size = 5) +
    ggtitle("Stage 10 – 12 Drosophila Embryonic Cell Type Labels") + 
    theme(text = element_text(size = 18))
  ggsave(file.path("cell_type_UMAP.png"), plot = p, width = 17, height = 10)
  
  p = DimPlot(wt_object, group.by = 'batch', label = FALSE, label.size = 5) +
    ggtitle("Stage 10 - 12: Batch") + 
    theme(text = element_text(size = 22), legend.text=element_text(size=22))
  ggsave(file.path("Batch_UMAP.png"), plot = p, width = 10, height = 8)
  
  p = DimPlot(wt_object, group.by = 'seurat_clusters', label = TRUE, label.size = 8) +
    ggtitle("Stage 10 - 12: Seurat Clusters") + 
    theme(text = element_text(size = 22), legend.text=element_text(size=22))
  ggsave(file.path("clusters_UMAP.png"), plot = p, width = 10, height = 8)
  
})

# plot out the proportion of each cell type in the stage 10-12 embryos 
withr::with_dir(TARGET_dir, {
  proportion_df = data.frame("cell_types" = names(table(wt_object@meta.data$manual_celltypes)), 
                             "number_cells" = as.vector(table(wt_object@meta.data$manual_celltypes)))
  proportion_df$cell_proportion = proportion_df$number_cells / sum(proportion_df$number_cells)
  
  p<-ggplot(data=proportion_df, aes(x=reorder(cell_types, cell_proportion), y=cell_proportion, fill = cell_types)) +
    scale_fill_discrete(name = "Cell Types") +
    geom_bar(stat="identity") + theme_bw() + coord_flip() + 
    ylab("Stage 10 - 12: total cell proportion") + 
    xlab("Cell Types") + 
    theme(text = element_text(size = 24), legend.position="none")
  ggsave(filename = file.path("cell_proportion_bar.png"), plot = p, width = 10, height = 10)
})


