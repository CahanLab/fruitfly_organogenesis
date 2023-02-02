library(ggplot2)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "match_early_clusters_to_late_clusters")
dir.create(TARGET_dir)

# look at spearman correlation between the two 
early_data = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))
late_data = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object3.rds"))

early_data@meta.data$celltype_cluster = paste0(early_data@meta.data$manual_celltypes, "_", early_data@meta.data$seurat_clusters)
late_data@meta.data$celltype_cluster = paste0(late_data@meta.data$manual_celltypes, "_", late_data@meta.data$seurat_clusters)

calc_spearman_correlation <- function(other_data = early_data, our_data = late_data) {
  other_meta = other_data@meta.data
  other_exp = other_data@assays$RNA@data
  
  our_meta = our_data@meta.data
  our_exp = our_data@assays$RNA@data
  
  combination_df = expand.grid(unique(our_meta$celltype_cluster), unique(other_meta$celltype_cluster))
  colnames(combination_df) = c('our_ct', 'other_ct')  
  combination_df$spearman_correlation = NA
  
  intersect_genes = intersect(rownames(other_exp), rownames(our_exp))
  other_exp = other_exp[intersect_genes, ]
  our_exp = our_exp[intersect_genes, ]
  
  our_avg_matrix = matrix(data = NA, 
                          nrow = nrow(our_exp), 
                          ncol = length(unique(combination_df$our_ct)))
  colnames(our_avg_matrix) = unique(combination_df$our_ct)
  rownames(our_avg_matrix) = rownames(our_exp)
  
  for(ct in colnames(our_avg_matrix)) {
    sub_our_met = our_meta[our_meta$celltype_cluster == ct, ]
    sub_our_exp = our_exp[, rownames(sub_our_met)]
    our_avg = apply(sub_our_exp, FUN = mean, MARGIN = 1)
    our_avg_matrix[, ct] = our_avg
  }
  
  other_avg_matrix = matrix(data = NA, 
                            nrow = nrow(other_exp), 
                            ncol = length(unique(combination_df$other_ct)))
  colnames(other_avg_matrix) = unique(combination_df$other_ct)
  rownames(other_avg_matrix) = rownames(other_exp)
  
  for(ct in colnames(other_avg_matrix)) {
    sub_other_met = other_meta[other_meta$celltype_cluster == ct, ]
    sub_other_exp = other_exp[, rownames(sub_other_met)]
    other_avg = apply(sub_other_exp, FUN = mean, MARGIN = 1)
    other_avg_matrix[, ct] = other_avg
  }
  
  for(temp_index in rownames(combination_df)) {
    print(temp_index)
    our_ct = combination_df[temp_index, 'our_ct']
    other_ct = combination_df[temp_index, 'other_ct']
    our_avg = our_avg_matrix[, our_ct]
    other_avg = other_avg_matrix[, other_ct]
    
    combination_df[temp_index, "spearman_correlation"] = cor(our_avg, other_avg, method = 'spearman')
  }
  
  return(combination_df)
}

combination_df = calc_spearman_correlation(other_data = early_data, our_data = late_data)
colnames(combination_df) = c("late_cell_types", "early_cell_types", "spearman_correlation")
saveRDS(combination_df, file = file.path(TARGET_dir, 'spearman_correlation_wt13_early_wt12.rds'))

combination_df = readRDS(file.path(TARGET_dir, 'spearman_correlation_wt13_early_wt12.rds'))
combination_df$late_cell_types = as.character(combination_df$late_cell_types)
combination_df$late_cell_types = factor(combination_df$late_cell_types, levels=sort(unique(combination_df$late_cell_types)))

combination_df$early_cell_types = as.character(combination_df$early_cell_types)
combination_df$early_cell_types = factor(combination_df$early_cell_types, levels=sort(unique(combination_df$early_cell_types)))

combination_df = readRDS(file.path(TARGET_dir, 'spearman_correlation_wt13_early_wt12.rds'))
match_df = data.frame("early_cts" = unique(combination_df$early_cell_types), 
                      "best_alignment" = NA)
for(ct in match_df$early_cts) {
  sub_comb_df = combination_df[combination_df$early_cell_types == ct, ]
  best_align = sub_comb_df[which.max(sub_comb_df$spearman_correlation), 'late_cell_types']
  match_df[match_df$early_cts == ct, 'best_alignment'] = as.vector(best_align)
}
write.csv(match_df, file = file.path(TARGET_dir, "early_match_late.csv"))

match_df = data.frame("late_cts" = unique(combination_df$late_cell_types), 
                      "best_alignment" = NA)
for(ct in match_df$late_cts) {
  sub_comb_df = combination_df[combination_df$late_cell_types == ct, ]
  best_align = sub_comb_df[which.max(sub_comb_df$spearman_correlation), 'early_cell_types']
  match_df[match_df$late_cts == ct, 'best_alignment'] = as.vector(best_align)
}
write.csv(match_df, file = file.path(TARGET_dir, "late_match_early.csv"))

library(viridis)
p = ggplot(combination_df, aes(late_cell_types, early_cell_types, fill= spearman_correlation)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  ggtitle("Cluster similarity between early and late data")
ggsave(filename = file.path(TARGET_dir, 'early_late_correlation.png'), plot = p, width = 14, height = 6)
