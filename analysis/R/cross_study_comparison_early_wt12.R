# perform comparisons between our data and data from Calderon et al, Seroka et al. 
TARGET_dir = file.path("results", ANALYSIS_VERSION, "cross_study_comparison_early_wt12")
dir.create(TARGET_dir)

##### harmonize the cell types #####
# load in our data 
our_object = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_early_wt12/manual_celltype_object1.rds'))

# load in Seroka data and only select stage 12
Seroka_object = readRDS("accessory_data/Doe_Drosophila_Embryo_Atlas/script/curated_embryo_Doe.rds")
Seroka_object = Seroka_object[, Seroka_object$dataset == 'stg12']

# load in Calderon data
Calderon_object = readRDS("accessory_data/continuum_drosophila_embryonic_development_RNA/processed_data/continuum_exploration/10-12_celltyped.rds")
Calderon_object@meta.data$cell_type = stringr::str_trim(Calderon_object@meta.data$cell_type)

# load in the harmonization table of embryo cell types 
label_map = read.csv("accessory_data/cross_study_harmonization/early_wt_stage10-12.csv")

#' @title harmonization of data
#' @description
#' harmonize the different cell types across datasets 
#' @param seurat_object the seurat object 
#' @param author the last name of the first author 
#' @param label_map dataframe that maps all the different cell types into a unified cell type 
#' @return seurat object with a harmonized cell type column in the meta data 
add_harmonize_ct <- function(seurat_object, author = 'Peng', label_map) {
  if(author == 'Peng') {
    col_id = 'Our.Data'
    seurat_ct_id = 'manual_celltypes'
  } else if(author == 'Calderon') {
    col_id = 'Calderon.et.al..2022'
    seurat_ct_id = 'cell_type'
  } else if(author == 'Seroka') {
    col_id = 'Seroka.et.al..2022'
    seurat_ct_id = 'cell_type'
  }
  
  sub_label_map = label_map[label_map[, col_id] != '', ]
  rownames(sub_label_map) = sub_label_map$Unified.Label
  seurat_object$harmonized_celltypes = NA
  
  for(h_ct in rownames(sub_label_map)) {
    orig_ct = sub_label_map[h_ct, col_id]
    orig_ct = stringr::str_split(orig_ct, ";")[[1]]
    orig_ct = orig_ct[orig_ct != '']
    seurat_object@meta.data[seurat_object@meta.data[, seurat_ct_id] %in% orig_ct, 'harmonized_celltypes'] = h_ct
  }
  return(seurat_object)
}

# harmonize our data 
our_object = add_harmonize_ct(our_object, author = 'Peng', label_map = label_map)

# harmonize seroka data
Seroka_object = add_harmonize_ct(Seroka_object, author = 'Seroka', label_map = label_map)

# harmonize calderon data 
Calderon_object = add_harmonize_ct(Calderon_object, author = 'Calderon', label_map = label_map)

##### calculate the proportion of cell types across different data ######
# calculate the various proportion of cell types in our data 
our_proportion_df = data.frame(cell_types = names(table(our_object@meta.data$harmonized_celltypes)), 
                               number = as.vector(table(our_object@meta.data$harmonized_celltypes)), 
                               data_type = 'our data')
our_proportion_df$proportion = our_proportion_df$number / sum(our_proportion_df$number)
zero_proportion_df = data.frame(cell_types = setdiff(unique(label_map$Unified.Label), unique(our_proportion_df$cell_types)), 
                                number = 0, 
                                data_type = 'our data', 
                                proportion = 0)
our_proportion_df = rbind(our_proportion_df, zero_proportion_df)

# calculate the various proportion of cell types in Seroka data 
Seroka_proportion_df = data.frame(cell_types = names(table(Seroka_object@meta.data$harmonized_celltypes)), 
                                  number = as.vector(table(Seroka_object@meta.data$harmonized_celltypes)), 
                                  data_type = 'Seroka et al, 2022')
Seroka_proportion_df$proportion = Seroka_proportion_df$number / sum(Seroka_proportion_df$number)
zero_proportion_df = data.frame(cell_types = setdiff(unique(label_map$Unified.Label), unique(Seroka_proportion_df$cell_types)), 
                                number = 0, 
                                data_type = 'Seroka et al, 2022', 
                                proportion = 0)
Seroka_proportion_df = rbind(Seroka_proportion_df, zero_proportion_df)

# calculate the various proportion of cell types in Calderon data 
Calderon_proportion_df = data.frame(cell_types = names(table(Calderon_object@meta.data$harmonized_celltypes)), 
                                  number = as.vector(table(Calderon_object@meta.data$harmonized_celltypes)), 
                                  data_type = 'Calderon et al, 2022')
Calderon_proportion_df$proportion = Calderon_proportion_df$number / sum(Calderon_proportion_df$number)
zero_proportion_df = data.frame(cell_types = setdiff(unique(label_map$Unified.Label), unique(Calderon_proportion_df$cell_types)), 
                                number = 0, 
                                data_type = 'Calderon et al, 2022', 
                                proportion = 0)
Calderon_proportion_df = rbind(Calderon_proportion_df, zero_proportion_df)

total_plot_df = rbind(our_proportion_df, Seroka_proportion_df, Calderon_proportion_df)
p<-ggplot(data=total_plot_df, aes(x=reorder(cell_types, proportion), y=proportion, fill = data_type)) +
  geom_bar(stat="identity", position=position_dodge()) + theme_bw() + coord_flip() + 
  ylab("Total Cell Proportion") + 
  xlab("Cell Types") + 
  ggtitle("Stage 10-12 Cell Type Proportions")

ggsave(filename = file.path(TARGET_dir, 'comparison_cell_proportion.png'), plot = p, width = 14, height = 6)
saveRDS(total_plot_df, file = file.path(TARGET_dir, "stage_10_12_proportion.rds"))

##### NOT USED -- compare spearman correlation between cell typing results across different data sets #####

#' @title calculate spearman correlation between cells types in our data and other data 
#' @param other_data the seurat object for other datasets 
#' @param our_data the seurat object for our dataset 
#' @param other_col the column of harmonized cell types for the other data 
#' @param our_col the column of harmonized cell types for our data 
calc_spearman_correlation <- function(other_data, our_data, other_col, our_col) {
  other_meta = other_data@meta.data
  other_exp = other_data@assays$RNA@data
  
  our_meta = our_data@meta.data
  our_exp = our_data@assays$RNA@data
  
  combination_df = expand.grid(unique(our_meta[, our_col]), unique(other_meta[, other_col]))
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
    sub_our_met = our_meta[our_meta[, our_col] == ct, ]
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
    sub_other_met = other_meta[other_meta[, other_col] == ct, ]
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

seroka_correlation = calc_spearman_correlation(Seroka_object, our_object, other_col = 'harmonized_celltypes', our_col =  'harmonized_celltypes')
saveRDS(seroka_correlation, file = file.path(TARGET_dir, "seroka_correlation.rds"))

seroka_correlation$our_ct = factor(seroka_correlation$our_ct, levels=sort(as.vector(unique(seroka_correlation$our_ct))))
seroka_correlation$other_ct = factor(seroka_correlation$other_ct, levels=sort(as.vector(unique(seroka_correlation$other_ct))))

p = ggplot(seroka_correlation, aes(our_ct, other_ct, fill= spearman_correlation)) + 
  geom_tile() +
  xlab("Our Cell Types") +
  ylab("Seroka et al's Cell Types") +
  scale_fill_viridis(option = "D", discrete=FALSE) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  ggtitle("Stage 10-12: Correlation similarity between cell types' centroids")
ggsave(filename = file.path(TARGET_dir, 'Seroka_cell_correlation.png'), plot = p, width = 14, height = 6)


Calderon_correlation = calc_spearman_correlation(Calderon_object, our_object, other_col = 'harmonized_celltypes', our_col =  'harmonized_celltypes')
saveRDS(Calderon_correlation, file = file.path(TARGET_dir, "calderon_correlation.rds"))

Calderon_correlation$our_ct = factor(Calderon_correlation$our_ct, levels=sort(as.vector(unique(Calderon_correlation$our_ct))))
Calderon_correlation$other_ct = factor(Calderon_correlation$other_ct, levels=sort(as.vector(unique(Calderon_correlation$other_ct))))

p = ggplot(Calderon_correlation, aes(our_ct, other_ct, fill= spearman_correlation)) + 
  geom_tile() +
  xlab("Our Cell Types") +
  ylab("Calderon et al's Cell Types") +
  scale_fill_viridis(option = "D", discrete=FALSE) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  ggtitle("Stage 10-12: Correlation similarity between cell types' centroids")
ggsave(filename = file.path(TARGET_dir, 'Calderon_cell_correlation.png'), plot = p, width = 14, height = 6)

##### compare agreement between our data's cell types and Seroka et al's cell types via pySCN #####
# write the training data (Seroka et al)
dir.create(file.path(TARGET_dir, 'Seroka_comparison_scn'))
withr::with_dir(
  file.path(TARGET_dir, 'Seroka_comparison_scn'), 
  { 
    raw_expTrain = Seroka_object@assays$RNA@counts
    raw_stTrain = Seroka_object@meta.data
    write(colnames(raw_expTrain), file = "raw_train_colnames.txt")
    write(rownames(raw_expTrain), file = "raw_train_rownames.txt")
    Matrix::writeMM(raw_expTrain, "raw_train_exp.txt")
    write.csv(raw_stTrain, file = 'raw_meta_tab.csv')
  }
)
# write the query data (our data)
withr::with_dir(
  file.path(TARGET_dir, 'Seroka_comparison_scn'), 
  {
    query_exp = our_object@assays$RNA@counts
    write(colnames(query_exp), file = file.path("raw_query_colnames.txt"))
    write(rownames(query_exp), file = file.path("raw_query_rownames.txt"))
    Matrix::writeMM(query_exp, file.path("raw_query_exp.txt"))
  }
)

# in the terminal run the python scripts in analysis/results/v18/cross_study_early_comparison_wt12/Seroka_comparison_scn
# load in the pySCN classification results
SCN_classification = read.csv(file.path(TARGET_dir, 'Seroka_comparison_scn', 'SCN_classification.csv'), row.names = 1)
SCN_classification = SCN_classification[rownames(our_object@meta.data), ]
our_object@meta.data$seroka_scn = SCN_classification$SCN_class

#' @title calculate the proportion of cells in our cell type classified as cell types from other studies 
#' @description
#' To assess the agreement between our cell typing results and those from other studies 
#' @param our_object our seruat object with pySCN classification results 
#' @param our_ct_col harmonized cell type column name 
#' @param other_ct_col pySCN classification column name
calc_class_proportion <- function(our_object, our_ct_col, other_ct_col) {
  our_meta = our_object@meta.data
  combination_df = expand.grid(unique(our_meta[, our_ct_col]), unique(our_meta[, other_ct_col]))
  colnames(combination_df) = c('our_ct', 'other_ct')  
  combination_df$class_proportion = NA
  for(temp_index in rownames(combination_df)) {
    our_ct = as.character(combination_df[temp_index, 'our_ct'])
    other_ct = as.character(combination_df[temp_index, 'other_ct'])
    
    sub_our_meta = our_meta[our_meta[, our_ct_col] == our_ct, ]
    other_proportion_list = table(sub_our_meta[, other_ct_col]) / nrow(sub_our_meta)
  
    if(other_ct %in% names(other_proportion_list) == FALSE) {
      combination_df[temp_index, "class_proportion"] = 0
    }
    else {
      combination_df[temp_index, "class_proportion"] = other_proportion_list[other_ct]
    }
  }
  return(combination_df)
}

seroka_proportion = calc_class_proportion(our_object, our_ct_col = 'harmonized_celltypes', 'seroka_scn')
seroka_proportion$other_ct = as.vector(seroka_proportion$other_ct)
seroka_proportion[seroka_proportion$other_ct == 'rand', 'other_ct'] = 'Unknown (SCN rand)'

# check if all the cell types are in there 
# if a cell type other than "Unknown" is not there, add it in and set the number to 0. This is mainly for plotting 
setdiff(c(unique(Seroka_object@meta.data$harmonized_celltypes), 'Unknown (SCN rand)'), unique(seroka_proportion$other_ct))
# [1] "Unknown"

seroka_proportion$our_ct = factor(seroka_proportion$our_ct, levels=sort(as.vector(unique(seroka_proportion$our_ct))))
seroka_proportion$other_ct = factor(seroka_proportion$other_ct, levels=sort(as.vector(unique(seroka_proportion$other_ct))))
saveRDS(seroka_proportion, file = file.path(TARGET_dir, "seroka_proportion.rds"))

p = ggplot(seroka_proportion, aes(our_ct, other_ct, fill= class_proportion)) + 
  geom_tile() +
  xlab("Our Cell Types") +
  ylab("Seroka et al's Cell Types") +
  scale_fill_viridis(option = "D", discrete=FALSE) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  ggtitle("Stage 10-12: SCN Classification Proportion (Seroka et al)")
ggsave(filename = file.path(TARGET_dir, 'Seroka_cell_proportion.png'), plot = p, width = 14, height = 6)


##### compare agreement between our data's cell types and Calderon et al's cell types via pySCN #####
# write the training data dir.create(file.path(TARGET_dir, 'Calderon_comparison_scn'))
withr::with_dir(
  file.path(TARGET_dir, 'Calderon_comparison_scn'), 
  { 
    raw_expTrain = Calderon_object@assays$RNA@counts
    raw_stTrain = Calderon_object@meta.data
    write(colnames(raw_expTrain), file = "raw_train_colnames.txt")
    write(rownames(raw_expTrain), file = "raw_train_rownames.txt")
    Matrix::writeMM(raw_expTrain, "raw_train_exp.txt")
    write.csv(raw_stTrain, file = 'raw_meta_tab.csv')
  }
)
# write the query data 
withr::with_dir(
  file.path(TARGET_dir, 'Calderon_comparison_scn'), 
  {
    query_exp = our_object@assays$RNA@counts
    write(colnames(query_exp), file = file.path("raw_query_colnames.txt"))
    write(rownames(query_exp), file = file.path("raw_query_rownames.txt"))
    Matrix::writeMM(query_exp, file.path("raw_query_exp.txt"))
  }
)

# in the terminal run the python scripts in analysis/results/v18/cross_study_comparison_early_wt12/Calderon_comparison_scn
SCN_classification = read.csv(file.path(TARGET_dir, 'Calderon_comparison_scn', 'SCN_classification.csv'), row.names = 1)
SCN_classification = SCN_classification[rownames(our_object@meta.data), ]
our_object@meta.data$calderon_scn = SCN_classification$SCN_class

calderon_proportion = calc_class_proportion(our_object, our_ct_col = 'harmonized_celltypes', 'calderon_scn')
calderon_proportion$other_ct = as.vector(calderon_proportion$other_ct)
calderon_proportion[calderon_proportion$other_ct == 'rand', 'other_ct'] = 'Unknown (SCN rand)'

# check if all the cell types are in there 
setdiff(c(unique(Calderon_object@meta.data$harmonized_celltypes), 'Unknown (SCN rand)'), unique(calderon_proportion$other_ct))
# [1] "Unknown" "Unknown (SCN rand)"
# we are missing "Unknown (SCN rand)" which means none of the cells got classified as that 
# to be consistent in our plotting, we will artifically add it in with class_proportion of 0 

# add in zero unknown (SCN rand)
zero_df = data.frame(our_ct = unique(calderon_proportion$our_ct),
                     other_ct = 'Unknown (SCN rand)', 
                     class_proportion = 0)
calderon_proportion = rbind(calderon_proportion, zero_df)
calderon_proportion$our_ct = factor(calderon_proportion$our_ct, levels=sort(as.vector(unique(calderon_proportion$our_ct))))
calderon_proportion$other_ct = factor(calderon_proportion$other_ct, levels=sort(as.vector(unique(calderon_proportion$other_ct))))
saveRDS(calderon_proportion, file = file.path(TARGET_dir, "calderon_proportion.rds"))

p = ggplot(calderon_proportion, aes(our_ct, other_ct, fill= class_proportion)) + 
  geom_tile() +
  xlab("Our Cell Types") +
  ylab("Calderon et al's Cell Types") +
  scale_fill_viridis(option = "D", discrete=FALSE) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  ggtitle("Stage 10-12: SCN Classification Proportion (Calderon et al)")
ggsave(filename = file.path(TARGET_dir, 'Calderon_cell_proportion.png'), plot = p, width = 14, height = 6)

##### run the pySCN classifier trained using our data's cell type label and applied to Seroka's dataset #####
# idea is to see if we can identify the rare cell types in our dataset in other datasets 
# write the training data (our data ) 
dir.create(file.path(TARGET_dir, 'Our_data_comparison_Seroka'))
withr::with_dir(
  file.path(TARGET_dir, 'Our_data_comparison_Seroka'), 
  {
    raw_train_exp = our_object@assays$RNA@counts
    train_st = our_object@meta.data

    write(colnames(raw_train_exp), file = file.path("raw_train_colnames.txt"))
    write(rownames(raw_train_exp), file = file.path("raw_train_rownames.txt"))
    Matrix::writeMM(raw_train_exp, file.path("raw_train_exp.txt"))
    
    write.csv(train_st, file = 'raw_train_st.csv')
  }
)

# write the query data
withr::with_dir(
  file.path(TARGET_dir, 'Our_data_comparison_Seroka'), 
  { 
    raw_expQuery = Seroka_object@assays$RNA@counts
    write(colnames(raw_expQuery), file = "raw_query_colnames.txt")
    write(rownames(raw_expQuery), file = "raw_query_rownames.txt")
    Matrix::writeMM(raw_expQuery, "raw_query_exp.txt")
  }
)

# in the terminal run the python scripts in analysis/results/v18/cross_study_comparison_early_wt12/Our_data_comparison_Seroka
SCN_classification = read.csv(file.path(TARGET_dir, 'Our_data_comparison_Seroka', 'SCN_classification.csv'), row.names = 1)
SCN_classification = SCN_classification[rownames(Seroka_object@meta.data), ]
Seroka_object@meta.data = cbind(Seroka_object@meta.data, SCN_classification)
saveRDS(Seroka_object, file = file.path(TARGET_dir, 'reverse_seroka_SCN_object.rds'))

##### run the pySCN classifier trained using our data's cell type label and applied to Calderon's dataset #####
# idea is to see if we can identify the rare cell types in our dataset in other datasets 

dir.create(file.path(TARGET_dir, 'Our_data_comparison_Calderon'))
#use the same classifier that was used to classify seroka
withr::with_dir(
  file.path(TARGET_dir, 'Our_data_comparison_Calderon'), 
  { 
    raw_expQuery = Calderon_object@assays$RNA@counts
    write(colnames(raw_expQuery), file = "raw_query_colnames.txt")
    write(rownames(raw_expQuery), file = "raw_query_rownames.txt")
    Matrix::writeMM(raw_expQuery, "raw_query_exp.txt")
  }
)

# in the terminal run the python scripts in analysis/results/v18/cross_study_comparison_early_wt12/Our_data_comparison_Calderon
SCN_classification = read.csv(file.path(TARGET_dir, 'Our_data_comparison_Calderon', 'SCN_classification.csv'), row.names = 1)
SCN_classification = SCN_classification[rownames(Calderon_object@meta.data), ]
Calderon_object@meta.data = cbind(Calderon_object@meta.data, SCN_classification)
saveRDS(Calderon_object, file = file.path(TARGET_dir, 'reverse_calderon_SCN_object.rds'))
