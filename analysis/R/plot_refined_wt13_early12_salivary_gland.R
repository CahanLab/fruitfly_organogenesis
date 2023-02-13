library(monocle3)
library(ggplot2)
library(RColorBrewer)
library(dbplyr)

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

UMAP_coord[UMAP_coord$clusters == 2, 'clusters'] = "Earlier Salivary Gland Cells"
UMAP_coord[UMAP_coord$clusters == 1, 'clusters'] = "Later Salivary Gland Cells"

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = pseudotime)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_viridis_c(option = "plasma") + 
  guides(fill=guide_legend(title="pseudo-time"))
ggsave(filename = file.path(TARGET_dir, "pseudotime.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = batch)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set1')
ggsave(filename = file.path(TARGET_dir, "batch.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = clusters)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set2')
ggsave(filename = file.path(TARGET_dir, "cluster.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=reorder(batch, pseudotime), y=pseudotime, fill = batch)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set1') + 
  ylab("pseudotime") + 
  xlab("batch")
ggsave(filename = file.path(TARGET_dir, "violin_pseudotime.png"), plot = p, width = 8, height = 6)

# plot out the GSEA results for early  
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "early_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[!is.na(GSEA_results$padj), ]
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]

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
  theme_bw()
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
  theme_bw()
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
index_list = vector()
for(gene in target_genes) { 
  index_list = c(index_list, which(tolower(gene) == tolower(rownames(norm_exp))))
}

norm_exp = norm_exp[index_list, ]
intersecting_genes = rownames(norm_exp)
norm_exp = monocle3::normalized_counts(cds)
norm_exp = as.matrix(norm_exp)
norm_exp = norm_exp[intersecting_genes, ]
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
png(filename = file.path(TARGET_dir, paste0(term, "_dynamic_gene_heatmap.png")), height = 2000, width = 1000, res = 200)
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


plot_df = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                     scaled_exp =   apply(scaled_exp, MARGIN = 2, FUN = mean), 
                     gene = 'vesicle transport related genes (average)')
plot_df$pseudotime = (plot_df$pseudotime - min(plot_df$pseudotime)) / max(plot_df$pseudotime)
p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
  xlab("pseudotime") + 
  ylab("average expression") +
  ggtitle(paste0('Average Gene Expression in ', term)) +
  geom_line(color = RColorBrewer::brewer.pal(n = 4, 'Set2')[1]) + theme_bw() 
ggsave(file.path(TARGET_dir, paste0(term, "_dynamic_gene_line_avg.png")), plot = p, width = 8, height = 5)

# plot out the Cytoplasmic translation
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "monocle3_no_batch_correct_object.rds"))
term = 'cytoplasmic translation (GO:0002181)'
GSEA_results = read.csv(file.path(TARGET_dir, "sig_late_GSEA_results.csv"), row.names = 1)

target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))

norm_exp = monocle3::normalized_counts(cds)
norm_exp = as.matrix(norm_exp)
index_list = vector()
for(gene in target_genes) { 
  index_list = c(index_list, which(tolower(gene) == tolower(rownames(norm_exp))))
}

norm_exp = norm_exp[index_list, ]
intersecting_genes = rownames(norm_exp)
norm_exp = monocle3::normalized_counts(cds)
norm_exp = as.matrix(norm_exp)
norm_exp = norm_exp[intersecting_genes, ]
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
png(filename = file.path(TARGET_dir, paste0(term, "_dynamic_gene_heatmap.png")), height = 2000, width = 1000, res = 200)
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


plot_df = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                     scaled_exp =   apply(scaled_exp, MARGIN = 2, FUN = mean), 
                     gene = 'translation related genes (average)')
plot_df$pseudotime = (plot_df$pseudotime - min(plot_df$pseudotime)) / max(plot_df$pseudotime)
p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
  xlab("pseudotime") + 
  ylab("average expression") +
  ggtitle(paste0('Average Gene Expression in ', term)) +
  geom_line(color = RColorBrewer::brewer.pal(n = 4, 'Set2')[2]) + theme_bw() 
ggsave(file.path(TARGET_dir, paste0(term, "_dynamic_gene_line_avg.png")), plot = p, width = 8, height = 5)


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

