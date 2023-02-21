library(Seurat)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'quality_metric_4_batches')
dir.create(TARGET_dir, recursive = TRUE)

early_object = readRDS(file.path("results/", ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))
late_object = readRDS(file.path("results/", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object4.rds"))

output_stats <- function(object) {
  cell_counts = ncol(object)
  meta_tab = object@meta.data
  
  median_reads = median(meta_tab$nCount_RNA)
  median_genes = median(meta_tab$nFeature_RNA)
  
  mean_reads = mean(meta_tab$nCount_RNA)
  mean_genes = mean(meta_tab$nFeature_RNA)
  
  return(c(cell_counts, median_reads, median_genes, mean_reads, mean_genes))
}

quality_mat = matrix(ncol = 5, nrow = 4, data = NA)
rownames(quality_mat) = c('Early WT 1', 'Early WT 2', 'Late WT 1', 'Late WT 2')
colnames(quality_mat) = c("cell number", 
                          "median UMI counts", 
                          "median gene counts", 
                          "mean UMI counts", 
                          'mean gene counts')

quality_mat['Early WT 1', ] = output_stats(early_object[, early_object$batch == 'rep_1'])
quality_mat['Early WT 2', ] = output_stats(early_object[, early_object$batch == 'rep_2'])
quality_mat['Late WT 1', ] = output_stats(late_object[, late_object$batch == 'rep_1'])
quality_mat['Late WT 2', ] = output_stats(late_object[, late_object$batch == 'rep_3'])

write.csv(quality_mat, file.path(TARGET_dir, 'quality_matrix.csv'))
