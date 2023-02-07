library(Seurat)
library(monocle3)
library(presto)
library(pheatmap)
library(ggplot2)
library(tradeSeq)

set.seed(123)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland")
dir.create(TARGET_dir)

wt_late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object3.rds"))

##########################
# TODO the below will be changed as time goes on 
wt_early_object = readRDS(file.path("results", ANALYSIS_VERSION, "early_wt12_integrated/BDGP_automated_annotation_object.rds"))
sub_wt_early_object = subset(wt_early_object, subset = Integrated_tentativeCellType == 'salivary gland body primordium')
sub_wt_early_object$experimental_condition = 'early'
########################

sub_wt_late_object = subset(wt_late_object, subset = manual_celltypes == "Salivary Gland")
sub_wt_late_object$experimental_condition = 'late'

combined_ct_object = merge(sub_wt_early_object, sub_wt_late_object)
object = Seurat::CreateSeuratObject(combined_ct_object@assays$RNA@counts, project = 'salivary_glands')
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
           genes=c("eyg"),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_toe.png'), width = 8, height = 6)

cds <- cluster_cells(cds, resolution = 1e-2)
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1)

marker_test_res <- monocle3::top_markers(cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)

write.csv(marker_test_res, file = file.path(TARGET_dir, "top_markers_monocle.csv"))

rank_sum_results = presto::wilcoxauc(normalized_counts(cds), cds@clusters$UMAP$clusters)
write.csv(rank_sum_results, file = file.path(TARGET_dir, "rank_sum_test.csv"))

plot_cells(cds, color_cells_by = "Phase", label_cell_groups = FALSE, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_phase.png'), width = 8, height = 6)

cds <- learn_graph(cds, use_partition = FALSE)

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

plot_cells(cds,
           genes=c("eyg", "CrebA", 'rib', 'wbl', 'toe', 'sage', 'fkh', 'sano', 'pip', 'RpS17', 'RpL24', 'RpL41', 'RpL13A'),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)

# this is where we use TradeSeq to test association 
expRaw = expression_matrix[apply(expression_matrix, MARGIN = 1, FUN = function(x){return(sum(x > 0)) }) > 10, ]
pt<-as.data.frame(pseudotime(cds))
colnames(pt)<-"pseudotime"
cw<-as.matrix(rep(1,nrow(pt)))
rownames(cw)<-rownames(pt)

ts<-tradeSeq::fitGAM(as.matrix(expRaw),pseudotime=as.matrix(pt),cellWeights=cw)
saveRDS(ts, file = file.path(TARGET_dir, "tradeseq_fitgam_results.rds"))
ATres<-tradeSeq::associationTest(ts)
saveRDS(ATres, file = file.path(TARGET_dir, 'raw_associationTest.rds'))

####################################################
# calculate the GO enrichment of the early vs late 
library(enrichR)
cds = readRDS(file.path(TARGET_dir, "monocle3_no_batch_correct_object.rds"))
enrichR::setEnrichrSite("FlyEnrichr")

rank_sum_test = read.csv(file.path(TARGET_dir, 'rank_sum_test.csv'), row.names = 1)
rank_sum_test = rank_sum_test[!is.na(rank_sum_test$padj), ]
rank_sum_test = rank_sum_test[rank_sum_test$padj < 0.05, ]
rank_sum_test = rank_sum_test[rank_sum_test$logFC > 0.1, ]

startRes_early = rank_sum_test[rank_sum_test$group == 2, ]
startRes_late = rank_sum_test[rank_sum_test$group == 1, ]

rownames(startRes_early) = startRes_early$feature
rownames(startRes_late) = startRes_late$feature

enrichment_results = enrichR::enrichr(
  genes = rownames(startRes_early), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018"
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
biological_analysis = biological_analysis[order(biological_analysis$Adjusted.P.value), ]
rownames(biological_analysis) = seq(1, nrow(biological_analysis))
molecular_analysis = enrichment_results$GO_Molecular_Function_2018
write.csv(biological_analysis, file = file.path(TARGET_dir, 'sig_GO_biological_early.csv'))

enrichment_results = enrichR::enrichr(
  genes = rownames(startRes_late), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018"
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
molecular_analysis = enrichment_results$GO_Molecular_Function_2018
write.csv(biological_analysis, file = file.path(TARGET_dir, 'sig_GO_biological_late.csv'))

# finish off with transient genes 
# start with transient genes
ATres = readRDS(file.path(TARGET_dir, 'raw_associationTest.rds'))
ATres = ATres[!is.na(ATres$pvalue), ]
ATres$adj_p = p.adjust(ATres$pvalue, method = 'fdr')
ATres = ATres[ATres$adj_p < 0.05, ]

rank_sum_test = read.csv(file.path(TARGET_dir, 'rank_sum_test.csv'), row.names = 1)
rank_sum_test = rank_sum_test[!is.na(rank_sum_test$padj), ]
rank_sum_test = rank_sum_test[rank_sum_test$padj < 0.05, ]
rank_sum_test = rank_sum_test[abs(rank_sum_test$logFC) > 0.1, ]

ATres = ATres[setdiff(rownames(ATres), rank_sum_test$feature), ]

enrichment_results = enrichR::enrichr(
  genes = rownames(ATres), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018"
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
molecular_analysis = enrichment_results$GO_Molecular_Function_2018
write.csv(biological_analysis, file = file.path(TARGET_dir, 'sig_GO_biological_transient.csv'))

#########################################################
library(fgsea)
pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')
rank_sum_test = read.csv(file.path(TARGET_dir, 'rank_sum_test.csv'), row.names = 1)
sub_rank_sum_test = rank_sum_test[rank_sum_test$group == 1, ]
ranks <- sub_rank_sum_test$logFC
names(ranks) <- sub_rank_sum_test$feature
fgseaRes <- fgsea(pathways = pathway_list, 
                  stats = ranks,
                  minSize=10,
                  maxSize=500)

fgseaRes = data.frame(fgseaRes) 
fgseaRes = apply(fgseaRes,2,as.character)
fgseaRes = as.data.frame(fgseaRes)
fgseaRes = fgseaRes[!is.na(fgseaRes$padj), ]
write.csv(fgseaRes, file = file.path(TARGET_dir, 'late_gsea_results.csv'))

sub_rank_sum_test = rank_sum_test[rank_sum_test$group == 2, ]
ranks <- sub_rank_sum_test$logFC
names(ranks) <- sub_rank_sum_test$feature
fgseaRes <- fgsea(pathways = pathway_list, 
                  stats = ranks,
                  minSize=10,
                  maxSize=500)

fgseaRes = data.frame(fgseaRes)
fgseaRes = apply(fgseaRes,2,as.character)
fgseaRes = as.data.frame(fgseaRes)
fgseaRes = fgseaRes[!is.na(fgseaRes$padj), ]
write.csv(fgseaRes, file = file.path(TARGET_dir, 'early_gsea_results.csv'))
