library(ggplot2)
library(RColorBrewer)
library(ggdendroplot)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'wt_late_all_marker_genes')
dir.create(TARGET_dir, recursive = TRUE)

object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object4.rds"))

selected_marker_genes = read.csv("accessory_data/marker_genes/wt13_marker_genes.csv")
selected_marker_genes = selected_marker_genes[order(selected_marker_genes$annotation), ]

p = DotPlot(
  object,
  assay = NULL,
  features = unique(selected_marker_genes$marker_genes),
  group.by = 'manual_celltypes', col.max = 3
) 
p = p + coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(file.path(TARGET_dir, 'genes_dotplot.png'), plot = p, width = 12, height = 10)
