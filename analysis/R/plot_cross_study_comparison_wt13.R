# plot out the cross study comparisons results for stage 13-16 wildtypes 
# Fig. 7c
# Supp Fig 15b, 33, 34

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'cross_study_comparison_wt13')
dir.create(TARGET_dir, recursive = TRUE)

calderon_SCN_matrix = readRDS(file.path('results', ANALYSIS_VERSION, "cross_study_comparison_wt13/calderon_proportion.rds"))
p = ggplot(calderon_SCN_matrix, aes(our_ct, other_ct, fill= class_proportion)) + 
  geom_tile() +
  xlab("Stage 13-16 Embryonic Cell Types") +
  ylab("SCN Cell Types from Calderon et al (stage 14-16)") +
  scale_fill_viridis(option = "D", discrete=FALSE, limits=c(0,1)) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  guides(fill=guide_legend(title="Percent Classified")) +
  theme(text = element_text(size = 22), plot.title.position = "plot") +
  ggtitle("Stage 13-16: SCN Classification Proportion using cell types from Calderon et al")
ggsave(filename = file.path(TARGET_dir, 'Calderon_SCN_proportion.png'), plot = p, width = 14, height = 10)

seroka_SCN_matrix = readRDS(file.path('results', ANALYSIS_VERSION, "cross_study_comparison_wt13/seroka_proportion.rds"))
p = ggplot(seroka_SCN_matrix, aes(our_ct, other_ct, fill= class_proportion)) + 
  geom_tile() +
  xlab("Stage 13-16 Embryonic Cell Types") +
  ylab("SCN Cell Types from Seroka et al (stage 14-16)") +
  scale_fill_viridis(option = "D", discrete=FALSE, limits=c(0,1)) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  guides(fill=guide_legend(title="Percent Classified")) +
  theme(text = element_text(size = 22), plot.title.position = "plot") +
  ggtitle("Stage 13-16: SCN Classification Proportion using cell types from Seroka et al") 
ggsave(filename = file.path(TARGET_dir, 'Seroka_SCN_proportion.png'), plot = p, width = 14, height = 10)

# this is to plot out cell type across different studies 
total_plot_df = readRDS(file.path('results', ANALYSIS_VERSION, "cross_study_comparison_wt13/stage_13_16_proportion.rds"))
total_plot_df[total_plot_df$data_type == 'our data', 'data_type'] = 'Stage 13-16 Embryonic Data'
total_plot_df[total_plot_df$data_type == 'Calderon et al, 2022', 'data_type'] = 'Calderon et al (stage 14-16)'
total_plot_df[total_plot_df$data_type == 'Seroka et al, 2022', 'data_type'] = 'Seroka et al (stage 14-16)'

p<-ggplot(data=total_plot_df, aes(x=reorder(cell_types, proportion), y=proportion, fill = data_type)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) + theme_bw() + coord_flip() + 
  ylab("Total Cell Proportion") + 
  xlab("Harmonized Cell Types") + 
  ggtitle("Stage 13-16 Cell Type Proportions") +
  ylim(c(0, 0.3)) +
  theme(text = element_text(size = 24), legend.title=element_blank(), plot.title.position = "plot")
ggsave(filename = file.path(TARGET_dir, 'comparison_cell_proportion.png'), plot = p, width = 10, height = 9)

# plot the reverse SCN results for Seroka data 
dir.create(file.path(TARGET_dir, 'reverse_SCN_seroka'), recursive = TRUE)
reverse_seroka_object = readRDS(file.path('results', ANALYSIS_VERSION, "cross_study_comparison_wt13/reverse_seroka_SCN_object.rds"))

withr::with_dir(file.path(TARGET_dir, 'reverse_SCN_seroka'), {
  start_index = grep("harmonized_celltypes", colnames(reverse_seroka_object@meta.data)) + 1
  
  for(scn_ct in unique(reverse_seroka_object@meta.data$SCN_class)) {
    reverse_seroka_object@meta.data$cur_ct = NA
    reverse_seroka_object@meta.data[reverse_seroka_object@meta.data$SCN_class == scn_ct, 'cur_ct'] = scn_ct
    reverse_seroka_object@meta.data[reverse_seroka_object@meta.data$SCN_class != scn_ct, 'cur_ct'] = 'Other'
    
    color_label = c(RColorBrewer::brewer.pal(n = 3, name = 'Set1')[1], '#c4c0c0')
    
    names(color_label) = c(scn_ct, 'Other')
    p = DimPlot(reverse_seroka_object, group.by = 'cur_ct') + 
        scale_colour_manual(values = color_label) + 
        #ggtitle(paste0("Stage 14-16 from Seroka et al, 2022: SCN classified ", scn_ct)) + 
        ggtitle("") +  
      theme(text = element_text(size = 18), plot.title.position = "plot", legend.text = element_text(size=18)) + 
      guides(color = guide_legend(override.aes = list(size = 10)))
    ggsave(paste(scn_ct, "_umap.png"), width = 10, height = 8)
    
  }
})

# plot the reverse SCN results for Calderon data
dir.create(file.path(TARGET_dir, 'reverse_SCN_calderon'), recursive = TRUE)
reverse_calderon_object = readRDS(file.path('results', ANALYSIS_VERSION, "cross_study_comparison_wt13/reverse_calderon_SCN_object.rds"))

withr::with_dir(file.path(TARGET_dir, 'reverse_SCN_calderon'), {

  for(scn_ct in unique(reverse_seroka_object@meta.data$SCN_class)) {
    reverse_calderon_object@meta.data$cur_ct = NA
    reverse_calderon_object@meta.data[reverse_calderon_object@meta.data$SCN_class == scn_ct, 'cur_ct'] = scn_ct
    reverse_calderon_object@meta.data[reverse_calderon_object@meta.data$SCN_class != scn_ct, 'cur_ct'] = 'Other'
    
    color_label = c(RColorBrewer::brewer.pal(n = 3, name = 'Set1')[1], '#c4c0c0')
    names(color_label) = c(scn_ct, 'Other')
    p = DimPlot(reverse_calderon_object, group.by = 'cur_ct') + 
      scale_colour_manual(values = color_label) + 
      #ggtitle(paste0("Stage 14-16 from Calderon et al, 2022: SCN classified ", scn_ct)) + 
      ggtitle("") +  
      theme(text = element_text(size = 18), plot.title.position = "plot", legend.text = element_text(size=18)) + 
      guides(color = guide_legend(override.aes = list(size = 10)))
    ggsave(paste(scn_ct, "_umap.png"), width = 10, height = 8)
    
  }
})


