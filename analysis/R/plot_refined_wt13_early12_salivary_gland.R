library(monocle3)
library(ggplot2)
library(RColorBrewer)
library(dbplyr)


plot_heatmap <- function(cds, target_genes, bandwidth = 3) {
  norm_exp = monocle3::normalized_counts(cds)
  norm_exp = as.matrix(norm_exp)
  norm_exp = norm_exp[c(target_genes), ]
  # this will change 
  #norm_exp = norm_exp[apply(norm_exp, MARGIN = 1, FUN = max) > 1, ]
  pt = monocle3::pseudotime(cds)
  pt = data.frame(pseudotime = pt)
  plot_df = cbind(pt, t(norm_exp[, rownames(pt)]))
  smoothed_df = data.frame()
  for(gene in colnames(plot_df)) {
    if(gene == 'pseudotime') {
      next
    }
    else {
      yy = ksmooth(plot_df[, 'pseudotime'], plot_df[, gene], kernel="normal", bandwidth = bandwidth, x.points=plot_df[, 'pseudotime'])
      if(nrow(smoothed_df) == 0) {
        smoothed_df = data.frame('pseudotime' = yy$x)
      }
      smoothed_df[, gene] = yy$y
    }
  }
  smoothed_df$pseudotime = NULL
  smoothed_df = t(smoothed_df)
  scaled_exp = t(scale(t(smoothed_df)))
  sorted_genes = names(sort(apply(scaled_exp, MARGIN = 1, FUN = which.max)))
  scaled_exp = scaled_exp[sorted_genes, ]
  return(scaled_exp)
}

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'refined_wt13_early12_salivary_gland')
dir.create(TARGET_dir, recursive = TRUE)

rank_sum = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "rank_sum_test.csv"), row.names = 1)
rank_sum = rank_sum[rank_sum$logFC > 0, ]
rank_sum[rank_sum$group == 2, 'group'] = "Earlier Salivary Gland Cells"
rank_sum[rank_sum$group == 1, 'group'] = "Later Salivary Gland Cells"

write.csv(rank_sum, file = file.path(TARGET_dir, 'DE_genes.csv'))
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "monocle3_no_batch_correct_object.rds"))

UMAP_coord = cds@int_colData$reducedDims$UMAP
colnames(UMAP_coord) = c("UMAP_1", "UMAP_2")
UMAP_coord = as.data.frame(UMAP_coord)
UMAP_coord$clusters = as.vector(monocle3::clusters(cds))
UMAP_coord$pseudotime = as.vector(monocle3::pseudotime(cds))
UMAP_coord$batch = as.vector(cds@colData$batch)

UMAP_coord[UMAP_coord$batch == 'early_rep_1', 'batch'] = 'Early rep 1'
UMAP_coord[UMAP_coord$batch == 'early_rep_2', 'batch'] = 'Early rep 2'
UMAP_coord[UMAP_coord$batch == 'late_rep_1', 'batch'] = 'Late rep 1'
UMAP_coord[UMAP_coord$batch == 'late_rep_3', 'batch'] = 'Late rep 2'

UMAP_coord[UMAP_coord$clusters == 2, 'clusters'] = "Early Salivary Gland Cells"
UMAP_coord[UMAP_coord$clusters == 1, 'clusters'] = "Late Salivary Gland Cells"

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = pseudotime)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_viridis_c(option = "plasma") + 
  guides(fill=guide_legend(title="pseudo-time")) + 
  theme(text = element_text(size = 18))
ggsave(filename = file.path(TARGET_dir, "pseudotime.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = batch)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set1') + 
  theme(text = element_text(size = 18))
ggsave(filename = file.path(TARGET_dir, "batch.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = clusters)) +
  geom_point() + 
  guides(color=guide_legend(title="")) +
  theme_minimal() + 
  scale_color_brewer(palette = 'Set2') + 
  theme(text = element_text(size = 18))
ggsave(filename = file.path(TARGET_dir, "celltypes.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=reorder(batch, pseudotime), y=pseudotime, fill = batch)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set1') + 
  ylab("pseudotime") + 
  xlab("batch") + 
  theme(text = element_text(size = 18))
ggsave(filename = file.path(TARGET_dir, "violin_pseudotime.png"), plot = p, width = 8, height = 6)

# plot out the GSEA results for early  
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "early_gsea_results.csv"))
GSEA_results = GSEA_results[!is.na(GSEA_results$padj), ]
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
GSEA_results = GSEA_results[, -1]

write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_early_GSEA_results.csv"))

# here are the interesting results that are not overlapping and are related to salivary gland development
focus_gsea = c("mRNA splicing, via spliceosome (GO:0000398)", 
               "regulation of transcription, DNA-templated (GO:0006355)",
               "mitotic cytokinesis (GO:0000281)", 
               "negative regulation of translation (GO:0017148)", 
               "salivary gland morphogenesis (GO:0007435)", 
               "dorsal closure (GO:0007391)", 
               "Golgi vesicle transport (GO:0048193)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[1]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("Genesets enriched in earlier salivary gland cells") +
  theme_bw() + 
  theme(text = element_text(size = 14), plot.title.position = "plot") 
ggsave(filename = file.path(TARGET_dir, "early_SG_GSEA_results.png"), plot = p, width = 8, height = 6)

# plot out the GSEA results for later  
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "late_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[!is.na(GSEA_results$padj), ]
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]

write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_late_GSEA_results.csv"))

# remove neuropeptide signaling pathway (GO:0007218) because the leading edge is too small 

#sub_GSEA_results = GSEA_results[GSEA_results$pathway != 'neuropeptide signaling pathway (GO:0007218)', ]
GSEA_results = GSEA_results[order(GSEA_results$padj), ]
sub_GSEA_results = GSEA_results[1:5, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[2]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("Genesets enriched in later salivary gland cells") +
  theme_bw() + 
  theme(text = element_text(size = 14), plot.title.position = "plot")
  
ggsave(filename = file.path(TARGET_dir, "later_SG_GSEA_results.png"), plot = p, width = 8, height = 6)


#############################
# plot out the Golgi Vesicle gene expression 
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "monocle3_no_batch_correct_object.rds"))
term = 'Golgi vesicle transport (GO:0048193)'
GSEA_results = read.csv(file.path(TARGET_dir, "sig_early_GSEA_results.csv"), row.names = 1)

target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))

norm_exp = monocle3::normalized_counts(cds)
norm_exp = as.matrix(norm_exp)
norm_exp = norm_exp[c('CrebA', target_genes), ]
# this will change 
#norm_exp = norm_exp[apply(norm_exp, MARGIN = 1, FUN = max) > 1, ]
pt = monocle3::pseudotime(cds)
pt = data.frame(pseudotime = pt)
plot_df = cbind(pt, t(norm_exp[, rownames(pt)]))
smoothed_df = data.frame()
for(gene in colnames(plot_df)) {
  if(gene == 'pseudotime') {
    next
  }
  else {
    yy = ksmooth(plot_df[, 'pseudotime'], plot_df[, gene], kernel="normal", bandwidth = 3, x.points=plot_df[, 'pseudotime'])
    if(nrow(smoothed_df) == 0) {
      smoothed_df = data.frame('pseudotime' = yy$x)
    }
    smoothed_df[, gene] = yy$y
  }
}
smoothed_df$pseudotime = NULL
smoothed_df = t(smoothed_df)
scaled_exp = t(scale(t(smoothed_df)))
sorted_genes = names(sort(apply(scaled_exp, MARGIN = 1, FUN = which.max)))
no_tf_scaled = scaled_exp[sorted_genes[sorted_genes != 'CrebA'], ]
png(filename = file.path(TARGET_dir, paste0(term, "_dynamic_gene_heatmap.png")), height = 2000, width = 1000, res = 200)
pheatmap(no_tf_scaled, cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

convert_line_plot <- function(scaled_exp) {
  plot_df = data.frame()
  for(gene in rownames(scaled_exp)) {
    temp_plot = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                           scaled_exp = scaled_exp[gene, ], 
                           gene = gene)
    plot_df = rbind(plot_df, temp_plot)
  }
  
  plot_df$pseudotime = plot_df$pseudotime / max(plot_df$pseudotime)
  return(plot_df)
}

plot_df = convert_line_plot(scaled_exp)
tf_plot_df = plot_df[plot_df$gene == 'CrebA', ]

scaled_exp = scaled_exp[rownames(scaled_exp) != 'CrebA', ]
plot_df = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                     scaled_exp =   apply(scaled_exp, MARGIN = 2, FUN = mean), 
                     gene = 'Golgi vesicle transport related genes (average)')
plot_df$pseudotime = (plot_df$pseudotime - min(plot_df$pseudotime)) / max(plot_df$pseudotime)
plot_df = rbind(plot_df, tf_plot_df)

no_tf_plot = plot_df[plot_df$gene != 'CrebA', ]
p<-ggplot(no_tf_plot, aes(x=pseudotime, y=scaled_exp, group=gene)) +
  xlab("pseudotime") + 
  ylab("average expression") +
  ggtitle(paste0('Average Gene Expression in ', term)) +
  geom_line(color = RColorBrewer::brewer.pal(n = 4, 'Set2')[1]) + theme_bw() 
ggsave(file.path(TARGET_dir, paste0(term, "_dynamic_gene_line_avg.png")), plot = p, width = 8, height = 5)

p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
  xlab("pseudotime") + 
  ylab("smoothed and scaled expression") +
  ggtitle(paste0('genes in ', term)) +
  geom_line(aes(color=gene)) + theme_bw() 
ggsave(file.path(TARGET_dir, paste0(term, "_CrebA_dynamic_gene_line_avg.png")), plot = p, width = 12, height = 8)


# plot out the Cytoplasmic translation
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "monocle3_no_batch_correct_object.rds"))
term = 'cytoplasmic translation (GO:0002181)'
GSEA_results = read.csv(file.path(TARGET_dir, "sig_late_GSEA_results.csv"), row.names = 1)

target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))

norm_exp = monocle3::normalized_counts(cds)
norm_exp = as.matrix(norm_exp)
norm_exp = norm_exp[c('rib', target_genes), ]
# this will change 
#norm_exp = norm_exp[apply(norm_exp, MARGIN = 1, FUN = max) > 1, ]
pt = monocle3::pseudotime(cds)
pt = data.frame(pseudotime = pt)
plot_df = cbind(pt, t(norm_exp[, rownames(pt)]))
smoothed_df = data.frame()
for(gene in colnames(plot_df)) {
  if(gene == 'pseudotime') {
    next
  }
  else {
    yy = ksmooth(plot_df[, 'pseudotime'], plot_df[, gene], kernel="normal", bandwidth = 3, x.points=plot_df[, 'pseudotime'])
    if(nrow(smoothed_df) == 0) {
      smoothed_df = data.frame('pseudotime' = yy$x)
    }
    smoothed_df[, gene] = yy$y
  }
}
smoothed_df$pseudotime = NULL
smoothed_df = t(smoothed_df)
scaled_exp = t(scale(t(smoothed_df)))
sorted_genes = names(sort(apply(scaled_exp, MARGIN = 1, FUN = which.max)))
no_tf_scaled = scaled_exp[sorted_genes[sorted_genes != 'rib'], ]
png(filename = file.path(TARGET_dir, paste0(term, "_dynamic_gene_heatmap.png")), height = 2000, width = 1000, res = 200)
pheatmap(no_tf_scaled, cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

convert_line_plot <- function(scaled_exp) {
  plot_df = data.frame()
  for(gene in rownames(scaled_exp)) {
    temp_plot = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                           scaled_exp = scaled_exp[gene, ], 
                           gene = gene)
    plot_df = rbind(plot_df, temp_plot)
  }
  
  plot_df$pseudotime = plot_df$pseudotime / max(plot_df$pseudotime)
  return(plot_df)
}

plot_df = convert_line_plot(scaled_exp)
tf_plot_df = plot_df[plot_df$gene == 'rib', ]

scaled_exp = scaled_exp[rownames(scaled_exp) != 'rib', ]
plot_df = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                     scaled_exp =   apply(scaled_exp, MARGIN = 2, FUN = mean), 
                     gene = 'translation related genes (average)')
plot_df$pseudotime = (plot_df$pseudotime - min(plot_df$pseudotime)) / max(plot_df$pseudotime)
plot_df = rbind(plot_df, tf_plot_df)

no_tf_plot = plot_df[plot_df$gene != 'rib', ]
p<-ggplot(no_tf_plot, aes(x=pseudotime, y=scaled_exp, group=gene)) +
  xlab("pseudotime") + 
  ylab("average expression") +
  ggtitle(paste0('Average Gene Expression in ', term)) +
  geom_line(color = RColorBrewer::brewer.pal(n = 4, 'Set2')[2]) + theme_bw() 
ggsave(file.path(TARGET_dir, paste0(term, "_dynamic_gene_line_avg.png")), plot = p, width = 8, height = 5)

p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
  xlab("pseudotime") + 
  ylab("smoothed and scaled expression") +
  ggtitle(paste0('genes in ', term)) +
  geom_line(aes(color=gene)) + theme_bw() 
ggsave(file.path(TARGET_dir, paste0(term, "_rib_dynamic_gene_line_avg.png")), plot = p, width = 10, height = 8)


##################################################
early_DE_genes = read.csv("results/v18/early_wt12_enrichment/Salivary Gland/markers_genes.csv", row.names = 1)
late_DE_genes = read.csv("results/v18/wt13_enrichment/Salivary Gland/markers_genes.csv", row.names = 1)

late_DE_genes = late_DE_genes[late_DE_genes$p_val_adj < 0.05 & late_DE_genes$avg_log2FC > 0, ]
early_DE_genes = early_DE_genes[early_DE_genes$p_val_adj < 0.05 & early_DE_genes$avg_log2FC > 0, ]

late_DE_genes$symbol = rownames(late_DE_genes)
early_DE_genes$symbol = rownames(early_DE_genes)

late_DE_genes$type = 'late'
early_DE_genes$type = 'early'
combined_DE_genes = rbind(late_DE_genes, early_DE_genes)
combined_DE_genes[combined_DE_genes$symbol %in% combined_DE_genes$symbol[duplicated(combined_DE_genes$symbol)], 'type'] = 
  paste0(combined_DE_genes[combined_DE_genes$symbol %in% combined_DE_genes$symbol[duplicated(combined_DE_genes$symbol)], 'type'], "_", "both")

sub_type_rank_sum = read.csv(file.path(TARGET_dir, 'DE_genes.csv'), row.names = 1)
sub_type_rank_sum = sub_type_rank_sum[sub_type_rank_sum$padj < 0.05, ]

early_sum_test = sub_type_rank_sum[sub_type_rank_sum$group == 'Earlier Salivary Gland Cells', ]
early_sum_test = early_sum_test[early_sum_test$feature %in% combined_DE_genes[combined_DE_genes$type == 'early', 'symbol'], ]

heatmap_df = plot_heatmap(cds, early_sum_test$feature, bandwidth = 3)

png(filename = file.path(TARGET_dir, paste0("early_SG_dynamic_gene_heatmap.png")), height = 2500, width = 1000, res = 200)
pheatmap::pheatmap(heatmap_df, cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

late_sum_test = sub_type_rank_sum[sub_type_rank_sum$group == 'Later Salivary Gland Cells', ]
late_sum_test = late_sum_test[late_sum_test$feature %in% combined_DE_genes[combined_DE_genes$type == 'late', 'symbol'], ]

heatmap_df = plot_heatmap(cds, late_sum_test$feature, bandwidth = 3)

png(filename = file.path(TARGET_dir, paste0("late_SG_dynamic_gene_heatmap.png")), height = 1500, width = 1000, res = 200)
pheatmap::pheatmap(heatmap_df, cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

# look at TFs 
TF_tab = read.csv("accessory_data/Drosophila_TFs/all_candidates.csv", sep = '\t')
TF_tab = TF_tab[TF_tab$verdict_DNA_BD != "NO", ]

early_TFs = intersect(TF_tab$symbol, early_sum_test$feature)
late_TFs = intersect(TF_tab$symbol, late_sum_test$feature)

middle_sum_test = combined_DE_genes[grep("both", combined_DE_genes$type), ]
middle_TF = intersect(TF_tab$symbol, unique(middle_sum_test$symbol))

# this is for the dotplot 
all_TFs = c(early_TFs, middle_TF, late_TFs)
cds@colData$cell_type = NA
cds@colData[clusters(cds) == 2, 'cell_type'] = "Early Salivary Gland Cells"
cds@colData[clusters(cds) == 1, 'cell_type'] = "Late Salivary Gland Cells"
  
meta_tab = cds@colData
norm_data = normalized_counts(cds)

diff_count = c()
for(TF in all_TFs) {
  early_norm = norm_data[, rownames(meta_tab[meta_tab$cell_type == 'Early Salivary Gland Cells', ])]
  late_norm = norm_data[, rownames(meta_tab[meta_tab$cell_type == 'Late Salivary Gland Cells', ])]
  
  percent_early = sum(early_norm[TF, ] > 0) / ncol(early_norm)  
  percent_late = sum(late_norm[TF, ] > 0) / ncol(late_norm) 
  
  diff = abs(percent_early - percent_late)
  diff_count = c(diff_count, percent_late)
}

names(diff_count) = all_TFs
sort(diff_count) 


p = plot_genes_by_group(cds, markers = names(sort(diff_count)), norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
  xlab("Cell Types") + coord_flip()

ggsave(filename = file.path(TARGET_dir, 'dynamic_TF.png'), plot = p, width = 10, height = 5)
zmeta_tab = cds@colData
norm_data = normalized_counts(cds)
meta_tab$sens = norm_data['sens', ]
meta_tab$fkh = norm_data['fkh', ]
meta_tab$CrebA = norm_data['CrebA', ]
meta_tab$sage = norm_data['sage', ]
meta_tab$toe = norm_data['toe', ]
meta_tab$eyg = norm_data['eyg', ]
meta_tab$trh = norm_data['trh', ]
meta_tab$rib = norm_data['rib', ]

meta_tab = as.data.frame(meta_tab)
p <- ggplot(meta_tab, aes(x=cell_type, y=rib)) + 
  geom_violin()

heatmap_df = plot_heatmap(cds, early_TFs, bandwidth = 3)
png(filename = file.path(TARGET_dir, paste0("early_TF_SG_dynamic_gene_heatmap.png")), height = 500, width = 1000, res = 200)
pheatmap::pheatmap(heatmap_df, cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

##################################################

# look at TFs specific 
early_rank_sum = read.csv("results/v18/early_wt12_enrichment/Salivary Gland/markers_genes.csv", row.names = 1)
late_rank_sum = read.csv("results/v18/wt13_enrichment/Salivary Gland/markers_genes.csv", row.names = 1)
early_rank_sum = early_rank_sum[early_rank_sum$p_val_adj < 0.05 & early_rank_sum$avg_log2FC > 0, ]
late_rank_sum = late_rank_sum[late_rank_sum$p_val_adj < 0.05 & late_rank_sum$avg_log2FC > 0, ]

TF_tab = read.csv("accessory_data/Drosophila_TFs/all_candidates.csv", sep = '\t')
TF_tab = TF_tab[TF_tab$verdict_DNA_BD != "NO", ]
all_genes = unique(c(rownames(late_rank_sum), rownames(early_rank_sum)))
i_TFs = intersect(TF_tab$symbol, all_genes) # TFs that are specific for salivary glands based on DE genes 

pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')
sg_TFs = intersect(pathway_list[["salivary gland development (GO:0007431)"]], TF_tab$symbol)
target_genes = intersect(i_TFs, sg_TFs) # TFs that are previously known to be 

early_object = readRDS("results/v18/manual_annotation_early_wt12/manual_celltype_object1.rds")

# plot out the violin plot for all the TFs 
dir.create(file.path(TARGET_dir, "early_TF"))
withr::with_dir("", {
  VlnPlot(early_object, features = 'SoxN', group.by = 'manual_celltypes', pt.size = 0)
})


###############################################
# this is plot out the dynamically expressed TFs in salivary gland development
ATres = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", 'raw_associationTest.rds'))
ATres = ATres[!is.na(ATres$pvalue), ]
ATres$adj_p = p.adjust(ATres$pvalue, method = 'fdr')
ATres = ATres[ATres$adj_p < 0.05, ]
pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')

TF_tab = read.csv("accessory_data/Drosophila_TFs/all_candidates.csv", sep = '\t')
TF_tab = TF_tab[TF_tab$verdict_DNA_BD != "NO", ]
i_TFs = intersect(TF_tab$symbol, rownames(ATres))

sg_TFs = intersect(pathway_list[["salivary gland development (GO:0007431)"]], TF_tab$symbol)
target_genes = intersect(i_TFs, sg_TFs)

norm_exp = monocle3::normalized_counts(cds)
norm_exp = as.matrix(norm_exp)
norm_exp = norm_exp[target_genes, ]
# this will change 
#norm_exp = norm_exp[apply(norm_exp, MARGIN = 1, FUN = max) > 1, ]
pt = monocle3::pseudotime(cds)
pt = data.frame(pseudotime = pt)
plot_df = cbind(pt, t(norm_exp[, rownames(pt)]))
smoothed_df = data.frame()
for(gene in colnames(plot_df)) {
  if(gene == 'pseudotime') {
    next
  }
  else {
    yy = ksmooth(plot_df[, 'pseudotime'], plot_df[, gene], kernel="normal", bandwidth = 3, x.points=plot_df[, 'pseudotime'])
    if(nrow(smoothed_df) == 0) {
      smoothed_df = data.frame('pseudotime' = yy$x)
    }
    smoothed_df[, gene] = yy$y
  }
}
smoothed_df$pseudotime = NULL
smoothed_df = t(smoothed_df)
scaled_exp = t(scale(t(smoothed_df)))
sorted_genes = names(sort(apply(scaled_exp, MARGIN = 1, FUN = which.max)))
scaled_exp = scaled_exp[sorted_genes, ]
png(filename = file.path(TARGET_dir, paste0("SG_TF_dynamic_gene_heatmap.png")), height = 2000, width = 1000, res = 200)
pheatmap(scaled_exp[sorted_genes, ], cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

convert_line_plot <- function(scaled_exp) {
  plot_df = data.frame()
  for(gene in rownames(scaled_exp)) {
    temp_plot = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                           scaled_exp = scaled_exp[gene, ], 
                           gene = gene)
    plot_df = rbind(plot_df, temp_plot)
  }
  
  plot_df$pseudotime = plot_df$pseudotime / max(plot_df$pseudotime)
  return(plot_df)
}


plot_df = convert_line_plot(scaled_exp)
plot_df$pseudotime = (plot_df$pseudotime - min(plot_df$pseudotime)) / max(plot_df$pseudotime)
rank_sum_test = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", 'rank_sum_test.csv'), row.names = 1)
rank_sum_test = rank_sum_test[rank_sum_test$padj < 0.05, ]
rank_sum_test = rank_sum_test[rank_sum_test$group == 1, ]
rank_sum_test = rank_sum_test[rank_sum_test$feature %in% plot_df$gene, ]
rank_sum_test = rank_sum_test[order(abs(rank_sum_test$logFC), decreasing = TRUE), ]

# interesting genes 
list_A = as.vector(rank_sum_test$feature[1:5])
sub_plot_df = plot_df[plot_df$gene %in% list_A, ]
p<-ggplot(sub_plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
  xlab("pseudotime") + 
  ylab("average expression") +
  geom_line(aes(color = gene)) + theme_bw() 
ggsave(file.path(TARGET_dir, "list_A_dynamic_gene_line_avg.png"), plot = p, width = 8, height = 5)

