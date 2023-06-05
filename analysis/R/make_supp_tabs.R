library(openxlsx)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'supplementary_tabs')
dir.create(TARGET_dir, recursive = TRUE)

##### make the marker genes with the corresponding cell types #####
marker_gene_list = list()
marker_gene_tab = read.csv("results/v18/wt13_integrated/marker_genes.csv")
marker_gene_tab$X = NULL
marker_gene_tab$cell_types = NA

seurat_object = readRDS("results/v18/manual_annotation_wt13/manual_celltype_object4.rds")
seurat_meta = seurat_object@meta.data
for(temp_cluster in unique(marker_gene_tab$cluster)) {
  marker_gene_tab[marker_gene_tab$cluster == temp_cluster, 'cell_types'] = unique(seurat_meta[seurat_meta$seurat_clusters == temp_cluster, "manual_celltypes"]) 
}

marker_gene_list[['stage 13-16']] = marker_gene_tab
# write in the stage 10-12 
marker_gene_tab = read.csv("results/v18/early_wt12_integrated/marker_genes.csv")
marker_gene_tab$X = NULL
marker_gene_tab$cell_types = NA

seurat_object = readRDS("results/v18/manual_annotation_early_wt12/manual_celltype_object1.rds")
seurat_meta = seurat_object@meta.data
for(temp_cluster in unique(marker_gene_tab$cluster)) {
  marker_gene_tab[marker_gene_tab$cluster == temp_cluster, 'cell_types'] = unique(seurat_meta[seurat_meta$seurat_clusters == temp_cluster, "manual_celltypes"]) 
}
marker_gene_list[['stage 10-12']] = marker_gene_tab

write.xlsx(marker_gene_list, file = file.path(TARGET_dir, "Table_2.xlsx"))

##### plot out the SG, trachea and GC #####
run_file_list = list()
run_file_list[['Stage13-16_SG']] = file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Salivary Gland/gsea_results_wt.csv")
run_file_list[['Stage13-16_Tr']] = file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Trachea/gsea_results_wt.csv")
run_file_list[['Stage13-16_GC']] = file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Germ Cells/gsea_results_wt.csv")
run_file_list[['Stage10-12_SG']] = file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Salivary Gland/gsea_results_wt.csv")
run_file_list[['Stage10-12_Tr']] = file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Trachea/gsea_results_wt.csv")
run_file_list[['Stage10-12_GC']] = file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Germ Cells/gsea_results_wt.csv")
table_list = list()
for(temp_name in names(run_file_list)) {
  gsea_results = read.csv(run_file_list[[temp_name]], row.names = 1)
  gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
  table_list[[temp_name]] = gsea_results
}
write.xlsx(table_list, file = file.path(TARGET_dir, "Table_3.xlsx"))

##### make Table 4 #####
run_file_list = list()
run_file_list[['Early_SG']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "early_gsea_results.csv")
run_file_list[['Late_SG']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "late_gsea_results.csv")
table_list = list()
for(temp_name in names(run_file_list)) {
  gsea_results = read.csv(run_file_list[[temp_name]])
  gsea_results$X = NULL
  gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
  table_list[[temp_name]] = gsea_results
}
write.xlsx(table_list, file = file.path(TARGET_dir, "Table_4.xlsx"))

##### make Table 5 #####
run_file_list = list()
run_file_list[['Tip_Tr']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Branching Trachea Cells_gsea_results.csv")
run_file_list[['Early_Tr']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Early Trachea Cells_gsea_results.csv")
run_file_list[['Interm_Tr']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Middle Trachea Cells_gsea_results.csv")
run_file_list[['Late_Tr']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Late Trachea Cells_gsea_results.csv")

table_list = list()
for(temp_name in names(run_file_list)) {
  gsea_results = read.csv(run_file_list[[temp_name]])
  gsea_results$X = NULL
  gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
  table_list[[temp_name]] = gsea_results
}
write.xlsx(table_list, file = file.path(TARGET_dir, "Table_5.xlsx"))

##### make Table 6 #####
run_file_list = list()
run_file_list[['Late_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "4_gsea_results.csv")
run_file_list[['Interm2_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "2_gsea_results.csv")
run_file_list[['Interm1_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "6_gsea_results.csv")
run_file_list[['Early_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "5_gsea_results.csv")
run_file_list[['Unknown1_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "1_gsea_results.csv")
run_file_list[['Unknown2_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "3_gsea_results.csv")

table_list = list()
for(temp_name in names(run_file_list)) {
  gsea_results = read.csv(run_file_list[[temp_name]])
  gsea_results$X = NULL
  gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
  table_list[[temp_name]] = gsea_results
}
write.xlsx(table_list, file = file.path(TARGET_dir, "Table_6.xlsx"))

##### make quality comparisons #####
output_statistics <- function(seurat_object) {
  return(c(median(seurat_object$nCount_RNA),
           median(seurat_object$nFeature_RNA), 
           mean(seurat_object$nCount_RNA),
           mean(seurat_object$nFeature_RNA)))
}

Seroka_object = readRDS("accessory_data/Doe_Drosophila_Embryo_Atlas/script/curated_embryo_Doe.rds")
Seroka_object = Seroka_object[, Seroka_object$dataset == 'stg12']
Calderon_object= readRDS("accessory_data/continuum_drosophila_embryonic_development_RNA/processed_data/continuum_exploration/10-12_celltyped.rds")
our_object = readRDS("results/v18/manual_annotation_early_wt12/manual_celltype_object1.rds")

stats_df = data.frame(row.names = c('median_nCount', 'median_nFeature', 'mean_nCount', 'mean_nFeature'))
stats_df$Seroka_stg12 = output_statistics(Seroka_object)
stats_df$Calderon_stg10_12 = output_statistics(Calderon_object)
stats_df$This_data_stg10_12 = output_statistics(our_object)

Seroka_object = readRDS("accessory_data/Doe_Drosophila_Embryo_Atlas/script/curated_embryo_Doe.rds")
Seroka_object = Seroka_object[, Seroka_object$dataset != 'stg12']
Calderon_object= readRDS("accessory_data/continuum_drosophila_embryonic_development_RNA/processed_data/continuum_exploration/14-16_celltyped.rds")
our_object = readRDS("results/v18/manual_annotation_wt13/manual_celltype_object4.rds")

stats_df$Seroka_stg14_16 = output_statistics(Seroka_object)
stats_df$Calderon_stg14_16 = output_statistics(Calderon_object)
stats_df$This_data_stg13_16 = output_statistics(our_object)

write.csv(stats_df, file = file.path(TARGET_dir, "Table_7.csv"))
