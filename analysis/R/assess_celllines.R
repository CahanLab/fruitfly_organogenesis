# correlate cell line bulk expression with cell types in embryos 
# Fig 7A
# Supp Fig 14
TARGET_dir = file.path("results", ANALYSIS_VERSION, 'identify_celllines')
dir.create(TARGET_dir)

##### load in object #####
object_wt13 = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object4.rds"))
object_early_wt12 = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))

###### get the marker genes across different cell types in both early and late embryos ######
all_marker_genes = c()
for(temp_path in list.dirs("results/v18/wt13_enrichment", recursive = FALSE)) {
  marker_genes = read.csv(file.path(temp_path, "markers_genes.csv"), row.names = 1)
  marker_genes = marker_genes[marker_genes$p_val_adj < 0.05, ]
  marker_genes = marker_genes[order(marker_genes$avg_log2FC, decreasing = TRUE), ]
  all_marker_genes = c(all_marker_genes, rownames(marker_genes)[1:60])
}
for(temp_path in list.dirs("results/v18/early_wt12_enrichment", recursive = FALSE)) {
  marker_genes = read.csv(file.path(temp_path, "markers_genes.csv"), row.names = 1)
  marker_genes = marker_genes[marker_genes$p_val_adj < 0.05, ]
  marker_genes = marker_genes[order(marker_genes$avg_log2FC, decreasing = TRUE), ]
  all_marker_genes = c(all_marker_genes, rownames(marker_genes)[1:60])
}
all_marker_genes = unique(all_marker_genes)

##### make the reference expression matrix #####
make_avg_tab <- function(object, all_marker_genes) {
  our_exp = object@assays$RNA@data
  our_exp = our_exp[all_marker_genes, ]
  
  our_meta = object@meta.data
  
  our_avg_matrix = matrix(data = NA, 
                          nrow = nrow(our_exp), 
                          ncol = length(unique(our_meta$manual_celltypes)))
  colnames(our_avg_matrix) = unique(our_meta$manual_celltypes)
  rownames(our_avg_matrix) = rownames(our_exp)
  
  for(ct in colnames(our_avg_matrix)) {
    sub_our_met = our_meta[our_meta$manual_celltypes == ct, ]
    sub_our_exp = our_exp[, rownames(sub_our_met)]
    our_avg = apply(sub_our_exp, FUN = mean, MARGIN = 1)
    our_avg_matrix[, ct] = our_avg
  }
  
  return(our_avg_matrix)
}

avg_tab_wt13 = make_avg_tab(object_wt13, all_marker_genes)
avg_tab_early_wt12 = make_avg_tab(object_early_wt12, all_marker_genes)

colnames(avg_tab_wt13) = paste0(colnames(avg_tab_wt13), "_stg13-16")
colnames(avg_tab_early_wt12) = paste0(colnames(avg_tab_early_wt12), "_stg10-12")

saveRDS(avg_tab_wt13, file = file.path(TARGET_dir, "avg_tab_wt13.rds"))
saveRDS(avg_tab_early_wt12, file = file.path(TARGET_dir, "avg_tab_early_wt12.rds"))

big_avg_tab = cbind(avg_tab_wt13, avg_tab_early_wt12)

##### load in the RPKM of cell lines and perform correlation between cell lines and reference matrix #####
cell_lines_df = read.csv("accessory_data/Drosophila_CellLines_FlyBase/gene_rpkm_matrix_fb_2023_02_modified.tsv", sep = '\t', row.names = 1)
cell_lines_df = cell_lines_df[all_marker_genes, ]

combination_df = expand.grid(unique(colnames(big_avg_tab)), unique(colnames(cell_lines_df)))
colnames(combination_df) = c('our_ct', 'cell_lines')  
combination_df$spearman_correlation = NA
cell_lines_df[is.na(cell_lines_df)] = 0
for(temp_index in rownames(combination_df)) {
  print(temp_index)
  our_ct = combination_df[temp_index, 'our_ct']
  other_ct = combination_df[temp_index, 'cell_lines']
  our_avg = big_avg_tab[, our_ct]
  other_avg = cell_lines_df[, other_ct]
  
  combination_df[temp_index, "spearman_correlation"] = cor(our_avg, other_avg, method = 'spearman')
}

write.csv(combination_df, file = file.path(TARGET_dir, 'spearman_correlation_cl_ct.csv'))

combination_df = read.csv(file.path(TARGET_dir, 'spearman_correlation_cl_ct.csv'), row.names = 1)
combination_df$stages = stringr::str_split_fixed(combination_df$our_ct,pattern = "_", n = 2)[, 2]
combination_df$cell_types = stringr::str_split_fixed(combination_df$our_ct,pattern = "_", n = 2)[, 1]

em_cls = c("1182.4H_", "GM2_", "Kc167_", "S1_", "S3_")
filter_correlation_df <- function(combination_df, em_cls) {
  filtered_df = data.frame()
  for(cell_line in em_cls) {
    temp_comb_df = combination_df[grep(cell_line, combination_df$cell_lines), ]
    temp_comb_df$cell_line = cell_line
    filtered_df = rbind(filtered_df, temp_comb_df)
  }
  return(filtered_df)
}

my_filtered_df = filter_correlation_df(combination_df, em_cls)

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

scale_correlation_df <- function(my_filtered_df) {
  my_filtered_df$scale_spearman_correlation = NA
  returned_df = data.frame()
  for(ct in unique(my_filtered_df$cell_line)) {
    temp_filtered_df = my_filtered_df[my_filtered_df$cell_line == ct, ]
    scaled_correlation = scale_values(temp_filtered_df[, 'spearman_correlation'])
    temp_filtered_df$scale_spearman_correlation = scaled_correlation
    returned_df = rbind(returned_df, temp_filtered_df)
  }  
  return(returned_df)
}

scaled_df = scale_correlation_df(my_filtered_df)
saveRDS(scaled_df, file = file.path(TARGET_dir, 'scaled_df.rds'))

scaled_df = readRDS(file.path(TARGET_dir, 'scaled_df.rds'))
scaled_df$scale_spearman_correlation = as.numeric(scaled_df$scale_spearman_correlation)
scaled_df$cell_line = stringr::str_remove_all(scaled_df$cell_line, "_")
p <- ggplot(data= scaled_df, mapping = aes_string(y = 'cell_line', x = 'cell_types')) +
  geom_tile(mapping = aes_string(fill = 'scale_spearman_correlation')) +
  guides(fill = guide_colourbar(title = 'Scaled Spearman Correlation')) +
  scale_fill_viridis_c(option = "plasma") + 
  labs(
    x = '',
    y = 'Embryonic Cell Lines'
  ) + 
  theme_classic()  + 
  facet_grid(
    cols = vars(stages),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) +
  ggtitle("Drosophila Embryonic Cell Lines")
ggsave(filename = file.path(TARGET_dir, "Cell_lines.png"), plot = p, width = 15, height = 4)

# get the top correlated category 
scaled_df$top_cat = 'Other'
for(ct in unique(scaled_df$cell_line)) {
  scaled_df[scaled_df$scale_spearman_correlation == 1, 'top_cat'] = "Top Cell Type"
}

p <- ggplot(data= scaled_df, aes(y = cell_line, x = cell_types, fill = top_cat)) +
  geom_tile(color = "black") +
  guides(fill = guide_legend(title = 'Top Correlated Cell Type')) +
  scale_fill_brewer(palette = 'Set2') +
  labs(
    x = '',
    y = 'Embryonic Cell Lines'
  ) + 
  theme_classic()  + 
  facet_grid(
    cols = vars(stages),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) +
  ggtitle("Drosophila Embryonic Cell Lines")
ggsave(filename = file.path(TARGET_dir, "Cell_lines_top.png"), plot = p, width = 15, height = 3)

##### plotting out the adult correlation #####
adult_tissue = combination_df[grepl("RNA.Seq_Profile_FlyAtlas2_", combination_df$cell_lines), ]
adult_tissue$tissue = stringr::str_remove_all(adult_tissue$cell_lines, "RNA.Seq_Profile_FlyAtlas2_")
adult_tissue$tissue = stringr::str_remove_all(adult_tissue$tissue, "_\\..*")

adult_tissue = scale_correlation_df(adult_tissue)
exclude_list = c("Head", "Crop", "Carcass", "Whole", 'Eye', 'AccessoryGland', 'Spermathecum')
for(bad_ct in exclude_list) {
  adult_tissue = adult_tissue[grepl(bad_ct, adult_tissue$tissue) == FALSE, ]
}


p <- ggplot(data= adult_tissue, mapping = aes_string(y = 'tissue', x = 'cell_types')) +
  geom_tile(mapping = aes_string(fill = 'scale_spearman_correlation')) +
  #scale.func(range = c(0, 100), limits = c(scale.min, scale.max)) +
  guides(fill = guide_colourbar(title = 'Scaled Spearman Correlation')) +
  scale_fill_viridis_c(option = "plasma") + 
  labs(
    x = '',
    y = 'Adult Tissues'
  ) + 
  theme_classic()  + 
  facet_grid(
    cols = vars(stages),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) +
  ggtitle("Drosophila Adult Tissues")
ggsave(filename = file.path(TARGET_dir, "adult_tissues.png"), plot = p, width = 15, height = 5)

adult_tissue$top_cat = 'Other'
adult_tissue[adult_tissue$scale_spearman_correlation == 1, 'top_cat'] = "Top Cell Type"

p <- ggplot(data= adult_tissue, aes(y = tissue, x = cell_types, fill = top_cat)) +
  geom_tile(color = "black") +
  #scale.func(range = c(0, 100), limits = c(scale.min, scale.max)) +
  guides(fill = guide_legend(title = 'Top Correlated Cell Type')) +
  scale_fill_brewer(palette = 'Set2') +
  labs(
    x = '',
    y = 'Embryonic Adult Tissue'
  ) + 
  theme_classic()  + 
  facet_grid(
    cols = vars(stages),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) +
  ggtitle("Drosophila Adult Tissues")
ggsave(filename = file.path(TARGET_dir, "Adult_tissue_top.png"), plot = p, width = 15, height = 5)

##### plotting out all the results #####
for(cell_line in unique(combination_df$cell_lines)) {
  temp_combo_df = combination_df[combination_df$cell_lines == cell_line, ]
  p<-ggplot(data=temp_combo_df, aes(x=reorder(our_ct, spearman_correlation), y=spearman_correlation)) +
    geom_bar(stat="identity") + coord_flip() + theme_bw() + 
    xlab("cell_types") + 
    ggtitle(cell_line)
  ggsave(file.path(TARGET_dir, paste0(stringr::str_replace_all(cell_line, "\\.", "_"), ".png")), height = 10, width = 8)
  
}
