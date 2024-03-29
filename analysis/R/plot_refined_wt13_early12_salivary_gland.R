# make plots for the analysis of salivary gland cells  
# Fig 3
# Supp Fig 3

##### load in data  ######
TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'refined_wt13_early12_salivary_gland')
dir.create(TARGET_dir, recursive = TRUE)

rank_sum = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "rank_sum_test.csv"), row.names = 1)
rank_sum = rank_sum[rank_sum$logFC > 0, ]
rank_sum[rank_sum$group == 2, 'group'] = "Earlier Salivary Gland Cells"
rank_sum[rank_sum$group == 1, 'group'] = "Later Salivary Gland Cells"

write.csv(rank_sum, file = file.path(TARGET_dir, 'DE_genes.csv'))
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "monocle3_no_batch_correct_object.rds"))

##### make umap and violin plots #####
UMAP_coord = cds@int_colData$reducedDims$UMAP
colnames(UMAP_coord) = c("UMAP_1", "UMAP_2")
UMAP_coord = as.data.frame(UMAP_coord)
UMAP_coord$clusters = as.vector(monocle3::clusters(cds))
UMAP_coord$pseudotime = as.vector(monocle3::pseudotime(cds))
UMAP_coord$batch = as.vector(cds@colData$batch)

UMAP_coord[UMAP_coord$batch == 'early_rep_1', 'batch'] = 'Stage 10-12 rep 1'
UMAP_coord[UMAP_coord$batch == 'early_rep_2', 'batch'] = 'Stage 10-12 rep 2'
UMAP_coord[UMAP_coord$batch == 'late_rep_1', 'batch'] = 'Stage 13-16 rep 1'
UMAP_coord[UMAP_coord$batch == 'late_rep_3', 'batch'] = 'Stage 13-16 rep 2'

UMAP_coord[UMAP_coord$clusters == 2, 'clusters'] = "Early Salivary Gland Cells"
UMAP_coord[UMAP_coord$clusters == 1, 'clusters'] = "Late Salivary Gland Cells"

# make umap plots with pseudotime
p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = pseudotime)) +
  geom_point(size = 3) + 
  theme_minimal() + 
  scale_color_viridis_c(option = "plasma") + 
  guides(fill=guide_legend(title="pseudo-time")) + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "pseudotime.png"), plot = p, width = 8, height = 6)

# make umap plots with batch information 
p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = batch)) +
  geom_point(size = 3) + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set1') + 
  theme(text = element_text(size = 24)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
ggsave(filename = file.path(TARGET_dir, "batch.png"), plot = p, width = 8, height = 6)

# make umap plot with cell types 
p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = clusters)) +
  geom_point(size = 3) + 
  guides(color=guide_legend(title="")) +
  theme_minimal() + 
  scale_color_brewer(palette = 'Set2') + 
  theme(text = element_text(size = 24)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
ggsave(filename = file.path(TARGET_dir, "celltypes.png"), plot = p, width = 10, height = 6)

# make violin plots with batch and pseudotime information
p = ggplot(UMAP_coord, aes(x=reorder(batch, pseudotime), y=pseudotime, fill = batch)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set1') + 
  ylab("Pseudotime") + 
  xlab("Batch") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_pseudotime.png"), plot = p, width = 8, height = 6)

###### plotting out GSEA results in bar plots ######
# plot out the GSEA results for early cells 
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "early_gsea_results.csv"))
GSEA_results = GSEA_results[!is.na(GSEA_results$padj), ]
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
GSEA_results = GSEA_results[, -1]

write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_early_GSEA_results.csv"))

# here are the interesting results that are not overlapping and are related to salivary gland development
focus_gsea = c("mRNA splicing, via spliceosome (GO:0000398)", 
               "regulation of transcription, DNA-templated (GO:0006355)",
               "negative regulation of translation (GO:0017148)", 
               "salivary gland morphogenesis (GO:0007435)", 
               "Golgi vesicle transport (GO:0048193)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot2::ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[1]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10(adj p-value)") + 
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 24), plot.title.position = "plot") 
ggsave(filename = file.path(TARGET_dir, "early_SG_GSEA_results.png"), plot = p, width = 10, height = 6)

# plot out the GSEA results for later cells  
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "late_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[!is.na(GSEA_results$padj), ]
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]

write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_late_GSEA_results.csv"))

GSEA_results = GSEA_results[order(GSEA_results$padj), ]
sub_GSEA_results = GSEA_results[1:5, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[2]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10(adj p-value)") + 
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 24), plot.title.position = "plot") 

ggsave(filename = file.path(TARGET_dir, "later_SG_GSEA_results.png"), plot = p, width = 10, height = 6)


##### making a function scale + smoothing expression profiles #####
#' @title make scaled and average expression pattern of genes across pseudotime 
#' @description
#' make scaled and smoothed expression profile of target genes across pseudotime. 
#' @param cds monocle3 object 
#' @param target_genes genes of interest 
#' @param bandwidth the bandwidth for ksmooth 
#' @return scaled expression profile
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

###### plotting out dynamic gene heatmap ######
# plot out the Golgi Vesicle gene expression 
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "monocle3_no_batch_correct_object.rds"))
term = 'Golgi vesicle transport (GO:0048193)'
GSEA_results = read.csv(file.path(TARGET_dir, "sig_early_GSEA_results.csv"), row.names = 1)

target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))

scaled_exp = plot_heatmap(cds, target_genes)
sorted_genes = names(sort(apply(scaled_exp, MARGIN = 1, FUN = which.max)))
no_tf_scaled = scaled_exp[sorted_genes[sorted_genes != 'CrebA'], ]
png(filename = file.path(TARGET_dir, paste0(term, "_dynamic_gene_heatmap.png")), height = 2000, width = 1000, res = 200)
pheatmap(no_tf_scaled, cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

# plot out the Cytoplasmic translation
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "monocle3_no_batch_correct_object.rds"))
term = 'cytoplasmic translation (GO:0002181)'
GSEA_results = read.csv(file.path(TARGET_dir, "sig_late_GSEA_results.csv"), row.names = 1)

target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))

norm_exp = monocle3::normalized_counts(cds)
norm_exp = as.matrix(norm_exp)
norm_exp = norm_exp[c('rib', target_genes), ]

scaled_exp = plot_heatmap(cds, target_genes)
sorted_genes = names(sort(apply(scaled_exp, MARGIN = 1, FUN = which.max)))
no_tf_scaled = scaled_exp[sorted_genes[sorted_genes != 'rib'], ]
png(filename = file.path(TARGET_dir, paste0(term, "_dynamic_gene_heatmap.png")), height = 2000, width = 1000, res = 200)
pheatmap(no_tf_scaled, cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

##### look at salivary gland specific genes #####
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "monocle3_no_batch_correct_object.rds"))
early_DE_genes = read.csv(file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Salivary Gland/markers_genes.csv"), row.names = 1)
late_DE_genes = read.csv(file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Salivary Gland/markers_genes.csv"), row.names = 1)

late_DE_genes = late_DE_genes[late_DE_genes$p_val_adj < 0.05 & late_DE_genes$avg_log2FC > 0, ]
early_DE_genes = early_DE_genes[early_DE_genes$p_val_adj < 0.05 & early_DE_genes$avg_log2FC > 0, ]

late_DE_genes$symbol = rownames(late_DE_genes)
early_DE_genes$symbol = rownames(early_DE_genes)

late_DE_genes$type = 'late'
early_DE_genes$type = 'early'
combined_DE_genes = rbind(late_DE_genes, early_DE_genes)
combined_DE_genes[combined_DE_genes$symbol %in% combined_DE_genes$symbol[duplicated(combined_DE_genes$symbol)], 'type'] = 
  paste0(combined_DE_genes[combined_DE_genes$symbol %in% combined_DE_genes$symbol[duplicated(combined_DE_genes$symbol)], 'type'], "_", "both")

combined_DE_genes = combined_DE_genes[combined_DE_genes$pct.1 > 0.15, ]
##### Plotting out the TFs #####
# look at TFs 
TF_tab = read.csv("accessory_data/Drosophila_TFs/all_candidates.csv", sep = '\t')
TF_tab = TF_tab[TF_tab$verdict_DNA_BD != "NO", ]

# this is for the dotplot 
all_TFs = unique(intersect(TF_tab$symbol, combined_DE_genes$symbol))
combined_DE_genes = combined_DE_genes[combined_DE_genes$symbol %in% all_TFs, ]
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

# dot plot
p = plot_genes_by_group(cds, markers = names(sort(diff_count)), norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
  xlab("Cell Types") + 
  ylab("Transcription Factors") +
  coord_flip() + 
  scale_x_discrete(limits = c('Late Salivary Gland Cells', 'Early Salivary Gland Cells')) + 
  theme(text = element_text(size = 24))

ggsave(filename = file.path(TARGET_dir, 'dynamic_TF.png'), plot = p, width = 15, height = 4.5)

