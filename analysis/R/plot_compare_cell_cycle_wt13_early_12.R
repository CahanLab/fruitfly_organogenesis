library(ggplot2)
library(RColorBrewer)
library(ggdendroplot)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'wt_late_cell_cycle')
dir.create(TARGET_dir, recursive = TRUE)

early_data = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))
late_data = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object3.rds"))

plot_df_late = data.frame(phase = names(table(late_data@meta.data$Phase) / sum(table(late_data@meta.data$Phase))),
                          proportion = as.vector(table(late_data@meta.data$Phase) / sum(table(late_data@meta.data$Phase))), 
                          data = 'late')

plot_df_early = data.frame(phase = names(table(early_data@meta.data$Phase) / sum(table(early_data@meta.data$Phase))),
                          proportion = as.vector(table(early_data@meta.data$Phase) / sum(table(early_data@meta.data$Phase))), 
                          data = 'early')

plot_df = rbind(plot_df_early, plot_df_late)

p = ggplot(data=plot_df, aes(x=phase, y=proportion, fill=data)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  guides(fill=guide_legend(title="Embryo Data")) +
  scale_fill_brewer(palette = "Set2") +
  ggtitle("All Cells in the Embryos") +
  ylim(c(0, 1))+
  theme_bw()

ggsave(filename = file.path(TARGET_dir, "all_cells_proportion.png"), plot = p, width = 8, height = 6)

intersecting_ct = intersect(unique(early_data@meta.data$manual_celltypes), unique(late_data@meta.data$manual_celltypes))

for(ct in intersecting_ct) {
  print(ct)
  sub_early_data = early_data[, early_data@meta.data[, 'manual_celltypes'] == ct]
  sub_late_data = late_data[, late_data@meta.data[, 'manual_celltypes'] == ct]
  
  plot_df_late = data.frame(phase = names(table(sub_late_data@meta.data$Phase) / sum(table(sub_late_data@meta.data$Phase))),
                            proportion = as.vector(table(sub_late_data@meta.data$Phase) / sum(table(sub_late_data@meta.data$Phase))), 
                            data = 'late')
  
  plot_df_early = data.frame(phase = names(table(sub_early_data@meta.data$Phase) / sum(table(sub_early_data@meta.data$Phase))),
                             proportion = as.vector(table(sub_early_data@meta.data$Phase) / sum(table(sub_early_data@meta.data$Phase))), 
                             data = 'early')
  
  plot_df = rbind(plot_df_early, plot_df_late)
  
  p = ggplot(data=plot_df, aes(x=phase, y=proportion, fill=data)) +
    geom_bar(stat="identity", position=position_dodge()) + 
    guides(fill=guide_legend(title="Embryo Data")) +
    scale_fill_brewer(palette = "Set2") +
    ggtitle(paste0(ct)) +
    ylim(c(0, 1)) +
    theme_bw()
  ggsave(filename = file.path(TARGET_dir, paste0(ct, "_cells_proportion.png")), plot = p, width = 8, height = 6)
  
}






