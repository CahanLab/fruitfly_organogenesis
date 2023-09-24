# plot out UMAPs of markers genes for SG, Tr, and SG 

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'UMAPS_marker_genes_early_3ct')
dir.create(TARGET_dir, recursive = TRUE)

object = readRDS(file.path('results', ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))

all_genes = c("CG14756", "Tpst", "Papss", "nur", "CG13159", "pip", "toe", "sage", "trh", 'Hs6st', 'CG3777', 'stumps', 'Gasp', 'Osi15', 'wisp', 'stai', 'CG4502', 'pgc', 'nos', 'Prosalpha5')
for(temp_gene in all_genes) {
  p = FeaturePlot(object, features = temp_gene) + 
    theme_void() + 
    ggtitle("") + 
    theme(legend.position = "none") 
  ggsave(filename = file.path(TARGET_dir, paste0(temp_gene, "_UMAP.png")), plot = p, width = 5, height = 5) 

}
