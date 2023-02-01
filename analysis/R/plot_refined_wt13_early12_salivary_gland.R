library(monocle3)
library(ggplot2)
library(RColorBrewer)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'refined_wt13_early12_salivary_gland')
dir.create(TARGET_dir, recursive = TRUE)

cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "monocle3_no_batch_correct_object.rds"))

UMAP_coord = cds@int_colData$reducedDims$UMAP
colnames(UMAP_coord) = c("UMAP_1", "UMAP_2")
UMAP_coord = as.data.frame(UMAP_coord)
UMAP_coord$clusters = as.vector(monocle3::clusters(cds))
UMAP_coord$pseudotime = as.vector(monocle3::pseudotime(cds))
UMAP_coord$batch = as.vector(cds@colData$batch)

UMAP_coord[UMAP_coord$batch == 'early_rep_1', 'batch'] = 'Early rep 1'
UMAP_coord[UMAP_coord$batch == 'early_rep_2', 'batch'] = 'Early rep 2'
UMAP_coord[UMAP_coord$batch == 'late_rep_1', 'batch'] = 'Late rep 1'
UMAP_coord[UMAP_coord$batch == 'late_rep_3', 'batch'] = 'Late rep 2'

UMAP_coord[UMAP_coord$clusters == 2, 'clusters'] = "Earlier Salivary Gland Cells"
UMAP_coord[UMAP_coord$clusters == 1, 'clusters'] = "Later Salivary Gland Cells"

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = pseudotime)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_viridis_c(option = "plasma") + 
  guides(fill=guide_legend(title="pseudo-time"))
ggsave(filename = file.path(TARGET_dir, "pseudotime.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = batch)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set2')
ggsave(filename = file.path(TARGET_dir, "batch.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = clusters)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set2')
ggsave(filename = file.path(TARGET_dir, "cluster.png"), plot = p, width = 8, height = 6)
