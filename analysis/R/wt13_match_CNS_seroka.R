library(Seurat)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt13_match_CNS_seroka")
dir.create(TARGET_dir)

seurat_data = readRDS("accessory_data/Doe_Drosophila_Embryo_Atlas/script/neuronal_subtype_Doe.rds")

our_data = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object3.rds"))
our_data = our_data[, our_data$manual_celltypes == 'CNS']

our_data$manual_celltypes = our_data$seurat_clusters
other_data = seurat_data
other_data$cell_type = seurat_data@active.ident

calc_spearman_correlation <- function(other_data, our_data) {
  other_meta = other_data@meta.data
  other_exp = other_data@assays$RNA@data
  
  our_meta = our_data@meta.data
  our_exp = our_data@assays$RNA@data
  
  combination_df = expand.grid(unique(our_meta$manual_celltypes), unique(other_meta$cell_type))
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
    sub_our_met = our_meta[our_meta$manual_celltypes == ct, ]
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
    sub_other_met = other_meta[other_meta$cell_type == ct, ]
    sub_other_exp = other_exp[, rownames(sub_other_met)]
    other_avg = apply(sub_other_exp, FUN = mean, MARGIN = 1)
    other_avg_matrix[, ct] = other_avg
  }
  
  for(temp_index in rownames(combination_df)) {
    print(temp_index)
    our_ct = combination_df[temp_index, 'our_ct']
    other_ct = combination_df[temp_index, 'other_ct']
    our_avg = our_avg_matrix[, as.vector(our_ct)]
    other_avg = other_avg_matrix[, other_ct]
    
    combination_df[temp_index, "spearman_correlation"] = cor(our_avg, other_avg, method = 'spearman')
  }
  
  return(combination_df)
}

combination_df = calc_spearman_correlation(other_data, our_data)
saveRDS(combination_df, file = file.path(TARGET_dir, 'spearman_correlation_wt13_Seroka.rds'))

library(viridis)
p = ggplot(combination_df, aes(our_ct, other_ct, fill= spearman_correlation)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  ggtitle("compare_neural_df")
ggsave(filename = file.path(TARGET_dir, 'seroka_cell_correlation.png'), plot = p, width = 14, height = 6)


p = DimPlot(our_data)
ggsave(filename = file.path(TARGET_dir, 'neural_UMAP.png'), plot = p, width = 14, height = 6)

Idents(our_data,  WhichCells(our_data, idents=33))<-'Enteric Neuron'
Idents(our_data,  WhichCells(our_data, idents=c(3, 11)))<-'newborn.N'
Idents(our_data,  WhichCells(our_data, idents=c(4, 15)))<-'old.N'

p = DimPlot(our_data)
ggsave(filename = file.path(TARGET_dir, 'neural_ct_UMAP.png'), plot = p, width = 14, height = 6)

p = FeaturePlot(our_data, features = c('mira', 'tap', 'Hey', 'nSyb', 'brp'))
ggsave(filename = file.path(TARGET_dir, 'neural_marker_genes_UMAP.png'), plot = p, width = 14, height = 10)

####################################
library(singleCellNet)

commonGenes = intersect(rownames(other_data@assays$RNA@counts), rownames(our_data@assays$RNA@counts))
exp_raw = other_data@assays$RNA@counts
st_raw = other_data@meta.data
exp_raw = exp_raw[commonGenes,]
set.seed(123)

stList = splitCommon(sampTab=st_raw, ncells=50, dLevel="cell_type")
stTrain = stList[[1]]
stTrain$barcode = rownames(stTrain)
expTrain = exp_raw[,rownames(stTrain)]
expTrain = as.matrix(expTrain)
system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 10, nRand = 70, nTrees = 1000, 
                                  nTopGenePairs = 25, dLevel = "cell_type", colName_samp = "barcode"))


query_exp = our_data@assays$RNA@counts
our_classification = scn_predict(class_info[['cnProc']], query_exp, nrand=12)

our_classification = our_classification[, rownames(our_data@meta.data)]
our_classification = t(our_classification)
our_data@meta.data = cbind(our_data@meta.data, our_classification)

p = FeaturePlot(our_data, c('NB', 'GMC', 'newborn.N', 'young.N', 'old.N'))
ggsave(filename = file.path(TARGET_dir, 'classification_scores_UMAP.png'), plot = p, width = 14, height = 10)
