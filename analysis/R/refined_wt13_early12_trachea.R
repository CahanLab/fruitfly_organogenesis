library(Seurat)
library(monocle3)
library(presto)
library(pheatmap)
library(ggplot2)
library(tradeSeq)

set.seed(123)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea")
dir.create(TARGET_dir)

wt_late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object3.rds"))

##########################
# TODO the below will be changed as time goes on 
wt_early_object = readRDS(file.path("results", ANALYSIS_VERSION, "early_wt12_integrated/BDGP_automated_annotation_object.rds"))
sub_wt_early_object = subset(wt_early_object, subset = Integrated_tentativeCellType == 'tracheal primordium')
sub_wt_early_object$experimental_condition = 'early'
########################

sub_wt_late_object = subset(wt_late_object, subset = manual_celltypes == "Trachea")
sub_wt_late_object$experimental_condition = 'late'

combined_ct_object = merge(sub_wt_early_object, sub_wt_late_object)
object = Seurat::CreateSeuratObject(combined_ct_object@assays$RNA@counts, project = 'Trachea')
object@meta.data$experimental_condition = combined_ct_object@meta.data$experimental_condition
object@meta.data$batch = paste0(combined_ct_object@meta.data$experimental_condition, "_", combined_ct_object@meta.data$batch)

object = Seurat::NormalizeData(object)
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

# start monocle3 
expression_matrix = object@assays$RNA@counts
cell_metadata = object@meta.data
gene_annotation = object@assays$RNA@meta.features
gene_annotation$gene_short_name = rownames(gene_annotation)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "batch", cell_size = 1, label_cell_groups = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_UMAP.png'), width = 8, height = 6)

plot_cells(cds,
           genes=c("Osi6"),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_toe.png'), width = 8, height = 6)

cds <- cluster_cells(cds, resolution = 5e-3)
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_cluster.png'), width = 8, height = 6)

marker_test_res <- monocle3::top_markers(cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)

write.csv(marker_test_res, file = file.path(TARGET_dir, "top_markers_monocle.csv"))
rank_sum_results = presto::wilcoxauc(normalized_counts(cds), cds@clusters$UMAP$clusters)
write.csv(rank_sum_results, file = file.path(TARGET_dir, "rank_sum_test.csv"))

plot_cells(cds, color_cells_by = "Phase", label_cell_groups = FALSE, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_phase.png'), width = 8, height = 6)

cds <- learn_graph(cds, use_partition = FALSE, learn_graph_control = list("minimal_branch_len" = 18))

plot_cells(cds,
           color_cells_by = "batch",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_pt.png'), width = 8, height = 6)
saveRDS(cds, file = file.path(TARGET_dir, "monocle3_no_batch_correct_object.rds"))

# this is where we use TradeSeq to test association 
expRaw = expression_matrix[apply(expression_matrix, MARGIN = 1, FUN = function(x){return(sum(x > 0) > 10) }), ]
pt<-as.data.frame(pseudotime(cds))
colnames(pt)<-"pseudotime"
cw<-as.matrix(rep(1,nrow(pt)))
rownames(cw)<-rownames(pt)
ts<-tradeSeq::fitGAM(as.matrix(expRaw),pseudotime=as.matrix(pt),cellWeights=cw, parallel = TRUE)
saveRDS(ts, file = file.path(TARGET_dir, "tradeseq_fitgam_results.rds"))

ATres<-associationTest(ts)
saveRDS(ATres, file = file.path(TARGET_dir, 'raw_associationTest.rds'))

plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1, cell_size = 1, show_trajectory_graph = FALSE)


cds@colData$cell_type = NA
cds@colData$cluster_id = monocle3::clusters(cds)
cds@colData[cds@colData$cluster_id %in% c(3, 7, 6), 'cell_type'] = 'Early Trachea Cells'
cds@colData[cds@colData$cluster_id %in% c(5, 4, 2), 'cell_type'] = 'Middle Trachea Cells'
cds@colData[cds@colData$cluster_id %in% c(8, 1), 'cell_type'] = 'Late Trachea Cells'
cds@colData[cds@colData$cluster_id %in% c(9), 'cell_type'] = 'Branching Trachea Cells'

saveRDS(cds, file = file.path(TARGET_dir, "monocle3_no_batch_correct_object.rds"))
##################################################
cds = readRDS(file.path(TARGET_dir, "monocle3_no_batch_correct_object.rds"))
rank_sum_results = presto::wilcoxauc(normalized_counts(cds), cds@colData$cell_type)
write.csv(rank_sum_results, file = file.path(TARGET_dir, "rank_sum_test.csv"))

library(fgsea) 
pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')
rank_sum_test = read.csv(file.path(TARGET_dir, 'rank_sum_test.csv'), row.names = 1)
rank_sum_test = rank_sum_test[rank_sum_test$pct_in > 10 | rank_sum_test$pct_out > 10, ]
for(ct in unique(rank_sum_test$group)) {
  sub_rank_sum_test = rank_sum_test[rank_sum_test$group == ct, ]
  ranks <- sub_rank_sum_test$logFC
  names(ranks) <- sub_rank_sum_test$feature
  fgseaRes <- fgsea(pathways = pathway_list, 
                    stats = ranks,
                    minSize=10,
                    maxSize=500)
  
  fgseaRes = data.frame(fgseaRes)
  fgseaRes = apply(fgseaRes,2,as.character)
  fgseaRes = as.data.frame(fgseaRes)
  fgseaRes$padj = as.numeric(fgseaRes$padj)
  fgseaRes = fgseaRes[!is.na(fgseaRes$padj), ]
  #fgseaRes = fgseaRes[fgseaRes$pval < 0.05, ]
  fgseaRes$NES = as.numeric(fgseaRes$NES)
  fgseaRes = fgseaRes[fgseaRes$NES > 0, ]
  write.csv(fgseaRes, file = file.path(TARGET_dir, paste0(ct, '_gsea_results.csv')))
}

# https://www.sciencedirect.com/science/article/pii/S1084952114000548#bib0090
# seems like Btl could be a diagnostic marker for the tip cells 
plot_cells(cds,
           genes=c('bnl', 'btl', 'esg', 'sty', 'Sec61alpha', 'Myo10A'), # 
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)


##################################################
# look at enrichr 
library(enrichR)
enrichR::setEnrichrSite("FlyEnrichr")
for(ct in unique(rank_sum_test$group)) {
  sub_rank_sum_test = rank_sum_test[rank_sum_test$group == ct, ]
  sub_rank_sum_test = sub_rank_sum_test[sub_rank_sum_test$logFC > 0.1, ]
  sub_rank_sum_test = sub_rank_sum_test[sub_rank_sum_test$padj < 0.05, ]
  
  enrichment_results = enrichR::enrichr(
    genes = sub_rank_sum_test$feature, 
    databases = c(
      "GO_Biological_Process_2018"
    )
  )
  
  biological_analysis = enrichment_results$GO_Biological_Process_2018
  
  write.csv(biological_analysis, file = file.path(TARGET_dir, paste0(ct, '_sig_GO_biological_early.csv')))
}



