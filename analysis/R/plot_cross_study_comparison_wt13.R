library(ggplot2)
library(RColorBrewer)
library(ggdendroplot)
library(viridis)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'cross_study_comparison_wt13')
dir.create(TARGET_dir, recursive = TRUE)

calderon_SCN_matrix = readRDS(file.path('results', ANALYSIS_VERSION, "cross_study_comparison_wt13/calderon_proportion.rds"))
p = ggplot(calderon_SCN_matrix, aes(our_ct, other_ct, fill= class_proportion)) + 
  geom_tile() +
  xlab("Stage 13-16 Embryonic Cell Types") +
  ylab("SCN Cell Types from Calderon et al (stage 14-16)") +
  scale_fill_viridis(option = "D", discrete=FALSE) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  guides(fill=guide_legend(title="Percent Classified")) +
  theme(text = element_text(size = 18)) +
  ggtitle("Stage 13-16: SCN Classification Proportion using cell types from Calderon et al")
ggsave(filename = file.path(TARGET_dir, 'Calderon_SCN_proportion.png'), plot = p, width = 14, height = 10)

seroka_SCN_matrix = readRDS(file.path('results', ANALYSIS_VERSION, "cross_study_comparison_wt13/seroka_proportion.rds"))
p = ggplot(seroka_SCN_matrix, aes(our_ct, other_ct, fill= class_proportion)) + 
  geom_tile() +
  xlab("Stage 13-16 Embryonic Cell Types") +
  ylab("SCN Cell Types from Seroka et al (stage 14-16)") +
  scale_fill_viridis(option = "D", discrete=FALSE) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  guides(fill=guide_legend(title="Percent Classified")) +
  theme(text = element_text(size = 18)) +
  ggtitle("Stage 13-16: SCN Classification Proportion using cell types from Seroka et al")
ggsave(filename = file.path(TARGET_dir, 'Seroka_SCN_proportion.png'), plot = p, width = 14, height = 10)

# this is to plot out cell type proportion 
total_plot_df = readRDS(file.path('results', ANALYSIS_VERSION, "cross_study_comparison_wt13/stage_13_16_proportion.rds"))
total_plot_df[total_plot_df$data_type == 'our data', 'data_type'] = 'Stage 13-16 Embryonic Data'
total_plot_df[total_plot_df$data_type == 'Calderon et al, 2022', 'data_type'] = 'Calderon et al (stage 14-16)'
total_plot_df[total_plot_df$data_type == 'Seroka et al, 2022', 'data_type'] = 'Seroka et al (stage 14-16)'

p<-ggplot(data=total_plot_df, aes(x=reorder(cell_types, proportion), y=proportion, fill = data_type)) +
  geom_bar(stat="identity", position=position_dodge()) + theme_bw() + coord_flip() + 
  ylab("Total Cell Proportion") + 
  xlab("Harmonized Cell Types") + 
  ggtitle("Stage 13-16 Cell Type Proportions") +
  theme(text = element_text(size = 18))

ggsave(filename = file.path(TARGET_dir, 'comparison_cell_proportion.png'), plot = p, width = 10, height = 14)




