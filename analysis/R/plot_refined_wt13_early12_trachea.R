library(monocle3)
library(ggplot2)
library(RColorBrewer)
library(dbplyr)


# Load in files -----------------------------------------------------------
# this section of code is there to load in all the appropriate files needed to make plots 
TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'refined_wt13_early12_trachea')
dir.create(TARGET_dir, recursive = TRUE)

rank_sum = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "rank_sum_test.csv"), row.names = 1)
rank_sum = rank_sum[rank_sum$logFC > 0, ]
write.csv(rank_sum, file = file.path(TARGET_dir, 'DE_genes.csv'))

cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "monocle3_no_batch_correct_object.rds"))

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

UMAP_coord$cell_type = cds@colData$cell_type
UMAP_coord[UMAP_coord$cell_type == 'Branching Trachea Cells', 'cell_type'] = 'Tracheal Tip Cells'
UMAP_coord[UMAP_coord$cell_type == 'Late Trachea Cells', 'cell_type'] = 'Late Tracheal Cells'
UMAP_coord[UMAP_coord$cell_type == 'Middle Trachea Cells', 'cell_type'] = 'Interm. Tracheal Cells'
UMAP_coord[UMAP_coord$cell_type == 'Early Trachea Cells', 'cell_type'] = 'Early Tracheal Cells'

# plotting umaps and violin plots -----------------------------------------
p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = pseudotime)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_viridis_c(option = "plasma") + 
  guides(fill=guide_legend(title="pseudo-time")) + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "pseudotime.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = batch)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set1') + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "batch.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = clusters)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set3') + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "cluster.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=reorder(batch, pseudotime), y=pseudotime, fill = batch)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set1') + 
  ylab("pseudotime") + 
  xlab("batch") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_pseudotime.png"), plot = p, width = 8, height = 6)


UMAP_coord$cell_type = factor(UMAP_coord$cell_type, levels = c('Early Tracheal Cells', 'Late Tracheal Cells', 'Interm. Tracheal Cells', 'Tracheal Tip Cells'))
p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = cell_type)) +
  guides(color=guide_legend(title="")) +
  geom_point() + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set2', breaks=c('Early Tracheal Cells', 'Interm. Tracheal Cells', 'Late Tracheal Cells', 'Tracheal Tip Cells')) + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "celltypes.png"), plot = p, width = 10, height = 6)

norm_exp = monocle3::normalized_counts(cds)

UMAP_coord$trh = norm_exp['trh', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -trh), y=trh, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2', breaks=c('Early Tracheal Cells', 'Interm. Tracheal Cells', 'Late Tracheal Cells', 'Tracheal Tip Cells')) + 
  ylab("trh normalized expression") + 
  xlab("cell type") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_trh.png"), plot = p, width = 8, height = 6)

UMAP_coord$Osi6 = norm_exp['Osi6', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -Osi6), y=Osi6, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  #geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2', breaks=c('Early Tracheal Cells', 'Interm. Tracheal Cells', 'Late Tracheal Cells', 'Tracheal Tip Cells')) + 
  ylab("Osi6 normalized expression") + 
  xlab("cell type") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_Osi6.png"), plot = p, width = 8, height = 6)

UMAP_coord$Osi17 = norm_exp['Osi17', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -Osi17), y=Osi17, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  #geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2', breaks=c('Early Tracheal Cells', 'Interm. Tracheal Cells', 'Late Tracheal Cells', 'Tracheal Tip Cells')) + 
  ylab("Osi17 normalized expression") + 
  xlab("cell type") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_Osi6.png"), plot = p, width = 8, height = 6)

UMAP_coord$btl = norm_exp['btl', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -btl), y=btl, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2', breaks=c('Early Tracheal Cells', 'Interm. Tracheal Cells', 'Late Tracheal Cells', 'Tracheal Tip Cells')) + 
  ylab("btl normalized expression") + 
  xlab("cell type") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_btl.png"), plot = p, width = 10, height = 6)

UMAP_coord$Mipp1 = norm_exp['Mipp1', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -Mipp1), y=Mipp1, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2', breaks=c('Early Tracheal Cells', 'Interm. Tracheal Cells', 'Late Tracheal Cells', 'Tracheal Tip Cells')) + 
  ylab("Mipp1 normalized expression") + 
  xlab("cell type") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_Mipp1.png"), plot = p, width = 10, height = 6)

UMAP_coord$bnl = norm_exp['bnl', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -bnl), y=bnl, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2', breaks=c('Early Tracheal Cells', 'Interm. Tracheal Cells', 'Late Tracheal Cells', 'Tracheal Tip Cells')) + 
  ylab("bnl normalized expression") + 
  xlab("cell type") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_bnl.png"), plot = p, width = 8, height = 6)

UMAP_coord$pnt = norm_exp['pnt', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -pnt), y=pnt, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2', breaks=c('Early Tracheal Cells', 'Interm. Tracheal Cells', 'Late Tracheal Cells', 'Tracheal Tip Cells')) + 
  ylab("pnt normalized expression") + 
  xlab("cell type") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_pnt.png"), plot = p, width = 8, height = 6)

UMAP_coord$sty = norm_exp['sty', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -sty), y=sty, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2', breaks=c('Early Tracheal Cells', 'Interm. Tracheal Cells', 'Late Tracheal Cells', 'Tracheal Tip Cells')) + 
  ylab("sty normalized expression") + 
  xlab("cell type") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(TARGET_dir, "violin_sty.png"), plot = p, width = 8, height = 6)

UMAP_coord$shg = norm_exp['shg', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -shg), y=shg, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2', breaks=c('Early Tracheal Cells', 'Interm. Tracheal Cells', 'Late Tracheal Cells', 'Tracheal Tip Cells')) + 
  ylab("shg normalized expression") + 
  xlab("cell type") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(filename = file.path(TARGET_dir, "violin_shg.png"), plot = p, width = 8, height = 6)


# Plot out GSEA results ---------------------------------------------------
# plot out the tracheal tip cells 
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Branching Trachea Cells_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_branching_GSEA_results.csv"))

# here are the interesting results that are not overlapping and are related to salivary gland development
focus_gsea = c("cell morphogenesis involved in differentiation (GO:0000904)", 
               'axon guidance (GO:0007411)',
               "cell-cell adhesion mediated by cadherin (GO:0044331)", 
               "positive regulation of intracellular signal transduction (GO:1902533)", 
               "apical junction assembly (GO:0043297)", 
               "establishment of epithelial cell apical/basal polarity (GO:0045198)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results[sub_GSEA_results$pathway == "calcium-dependent cell-cell adhesion via plasma membrane cell adhesion molecules (GO:0016339)", 'pathway'] = 'calcium-dependent cell-cell adhesion \n via plasma membrane cell adhesion molecules (GO:0016339)'
sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[4]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "Tip_Cells_GSEA_results.png"), plot = p, width = 10.4, height = 5)

# plot out all the genes associated with the terms 
for(term in sub_GSEA_results$pathway) { 
  target_genes = sub_GSEA_results[sub_GSEA_results$pathway == term, 'leadingEdge']
  target_genes = eval(parse(text = target_genes))
  
  # dot plot 
  p = plot_genes_by_group(cds, markers = target_genes, norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
    xlab("Cell Types") + coord_flip()
  ggsave(filename = file.path(TARGET_dir, paste0(term, '_genes.png')), plot = p, width = 15, height = 5)
}

# plot out the tracheal early cells 
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Early Trachea Cells_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_early_GSEA_results.csv"))

# here are the interesting results that are not overlapping and are related to salivary gland development
focus_gsea = c("mRNA splicing, via spliceosome (GO:0000398)",
               "regulation of gene expression (GO:0010468)", 
               "regulation of mitotic cell cycle (GO:0007346)", 
               "tracheal outgrowth, open tracheal system (GO:0007426)", 
               "Notch signaling pathway (GO:0007219)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[1]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 24))

ggsave(filename = file.path(TARGET_dir, "Early_Cells_GSEA_results.png"), plot = p, width = 10, height = 4)

for(term in sub_GSEA_results$pathway) { 
  target_genes = sub_GSEA_results[sub_GSEA_results$pathway == term, 'leadingEdge']
  target_genes = eval(parse(text = target_genes))
  
  # dot plot 
  p = plot_genes_by_group(cds, markers = target_genes, norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
    xlab("Cell Types") + coord_flip()
  ggsave(filename = file.path(TARGET_dir, paste0(term, '_genes.png')), plot = p, width = 15, height = 5)
}

# plot out the tracheal middle cells 
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Middle Trachea Cells_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_middle_GSEA_results.csv"))

# here are the interesting results that are not overlapping and are related to salivary gland development
focus_gsea = c("regulation of tube size, open tracheal system (GO:0035151)", 
               "septate junction assembly (GO:0019991)", 
               "regulation of translation (GO:0006417)", 
               "Golgi vesicle transport (GO:0048193)", 
               "chitin-based cuticle development (GO:0040003)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[3]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 24), plot.margin = margin(0,0,0,0.1, "cm"))

ggsave(filename = file.path(TARGET_dir, "Middle_Cells_GSEA_results.png"), plot = p, width = 10, height = 4)

for(term in sub_GSEA_results$pathway) { 
  target_genes = sub_GSEA_results[sub_GSEA_results$pathway == term, 'leadingEdge']
  target_genes = eval(parse(text = target_genes))
  
  # dot plot 
  p = plot_genes_by_group(cds, markers = target_genes, norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
    xlab("Cell Types") + coord_flip()
  ggsave(filename = file.path(TARGET_dir, paste0(term, '_genes.png')), plot = p, width = 15, height = 5)
}

# plot out the tracheal Late cells 
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Late Trachea Cells_gsea_results.csv"), row.names = 1)
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
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw() +   
  theme(text = element_text(size = 24))

ggsave(filename = file.path(TARGET_dir, "Late_Cells_GSEA_results.png"), plot = p, width = 8.8, height = 4)

for(term in sub_GSEA_results$pathway) { 
  target_genes = sub_GSEA_results[sub_GSEA_results$pathway == term, 'leadingEdge']
  target_genes = eval(parse(text = target_genes))
  
  # dot plot 
  p = plot_genes_by_group(cds, markers = target_genes, norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
    xlab("Cell Types") + coord_flip()
  ggsave(filename = file.path(TARGET_dir, paste0(term, '_genes.png')), plot = p, width = 15, height = 5)
}


# calculate cross correlation  --------------------------------------------
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

average_geneset <- function(scaled_exp) {
  plot_df = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                       scaled_exp =   apply(scaled_exp, MARGIN = 2, FUN = mean))
  plot_df$pseudotime = (plot_df$pseudotime - min(plot_df$pseudotime)) / max(plot_df$pseudotime)
  return(plot_df)
}
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "monocle3_no_batch_correct_object.rds"))
cds = cds[, monocle3::clusters(cds) != 9] # remove the branch from the main trajectory

GSEA_results = read.csv(file.path(TARGET_dir, "sig_middle_GSEA_results.csv"), row.names = 1)

term = 'regulation of tube size, open tracheal system (GO:0035151)'
target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
norm_exp_1 = plot_heatmap(cds, target_genes)
plot_norm_1 = average_geneset(norm_exp_1)
plot_norm_1$geneset = 'regulation of tube size (average expression)'

term = 'Golgi vesicle transport (GO:0048193)'
target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
norm_exp_2 = plot_heatmap(cds, target_genes)
plot_norm_2 = average_geneset(norm_exp_2)
plot_norm_2$geneset = 'Golgi vesicle transport (average expression)'

x = ccf(plot_norm_2$scaled_exp, plot_norm_1$scaled_exp)
correlation = x$acf[, 1, 1]
names(correlation) = x$lag[, 1, 1]
max_correlation = correlation[which.max(correlation)]

plot_df = rbind(plot_norm_1, plot_norm_2)
plot_df$geneset = stringr::str_replace_all(plot_df$geneset, "\\(", "\n\\(\\")

p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=geneset, color = geneset)) +
  xlab("pseudotime") + 
  ylab("average expression") +
  ggtitle(paste0('Golg vesicle transport vs regulation of tube size')) +
  geom_line() + 
  theme_bw() + 
  scale_color_brewer(palette='Set2') +
  annotate("text", x=0.8, y=1, label= paste0("cross-correlation: ", round(as.numeric(max_correlation), digits = 3),  "\n lag: ", round(as.numeric(names(max_correlation)) / nrow(plot_norm_1), digits = 3))) + 
  guides(color=guide_legend(title="")) + 
  theme(text = element_text(size = 18))
ggsave(file.path(TARGET_dir, paste0(term, "golgi_vs_tube_size.png")), plot = p, width = 8, height = 5)

# look at correlation between chitin and tube size
term = 'chitin-based cuticle development (GO:0040003)'
target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
norm_exp_2 = plot_heatmap(cds, target_genes)
plot_norm_2 = average_geneset(norm_exp_2)
plot_norm_2$geneset = 'chitin-based cuticle (average expression)'

x = ccf(plot_norm_2$scaled_exp, plot_norm_1$scaled_exp)
correlation = x$acf[, 1, 1]
names(correlation) = x$lag[, 1, 1]
max_correlation = correlation[which.max(correlation)]

plot_df = rbind(plot_norm_1, plot_norm_2)
plot_df$geneset = stringr::str_replace_all(plot_df$geneset, "\\(", "\n\\(\\")

p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=geneset, color = geneset)) +
  xlab("pseudotime") + 
  ylab("average expression") +
  ggtitle(paste0('chitin-based cuticle vs regulation of tube size')) +
  geom_line() + 
  theme_bw() + 
  scale_color_brewer(palette='Set2') +
  annotate("text", x=0.8, y=1, label= paste0("cross-correlation: ", round(as.numeric(max_correlation), digits = 3),  "\n lag: ", round(as.numeric(names(max_correlation)) / nrow(plot_norm_1), digits = 3))) + 
  guides(color=guide_legend(title="")) + 
  theme(text = element_text(size = 18))
ggsave(file.path(TARGET_dir, paste0(term, "chitin_vs_tube_size.png")), plot = p, width = 8, height = 5)

# plot out all three dynamics 
term = 'regulation of tube size, open tracheal system (GO:0035151)'
target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
norm_exp_1 = plot_heatmap(cds, target_genes)
plot_norm_1 = average_geneset(norm_exp_1)
plot_norm_1$geneset = 'regulation of tube size (average expression)'

term = 'Golgi vesicle transport (GO:0048193)'
target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
norm_exp_2 = plot_heatmap(cds, target_genes)
plot_norm_2 = average_geneset(norm_exp_2)
plot_norm_2$geneset = 'Golgi vesicle transport (average expression)'

term = 'chitin-based cuticle development (GO:0040003)'
target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
norm_exp_3 = plot_heatmap(cds, target_genes)
plot_norm_3 = average_geneset(norm_exp_3)
plot_norm_3$geneset = 'chitin-based cuticle development (average expression)'

x = ccf(plot_norm_2$scaled_exp, plot_norm_1$scaled_exp)
correlation = x$acf[, 1, 1]
names(correlation) = x$lag[, 1, 1]
max_correlation = correlation[which.max(correlation)]

plot_df = rbind(plot_norm_1, plot_norm_2, plot_norm_3)
plot_df$geneset = stringr::str_replace_all(plot_df$geneset, "\\(", "\n\\(\\")

p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=geneset, color = geneset)) +
  xlab("pseudotime") + 
  ylab("average expression") +
  ggtitle("") +
  geom_line(linewidth = 1.2) + 
  theme_bw() + 
  scale_color_brewer(palette='Set1') +
  guides(color=guide_legend(title="")) + 
  theme(text = element_text(size = 24), legend.position="bottom", plot.margin = margin(t = 0, r = 2, b = 0, l = 0, "cm"))
ggsave(file.path(TARGET_dir, paste0("golig_chitin_size.png")), plot = p, width = 11, height = 9)

# plot out the dynamic gene heatmap ---------------------------------------
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "monocle3_no_batch_correct_object.rds"))
cds = cds[, monocle3::clusters(cds) != 9] # remove the branch from the main trajectory

GSEA_results = read.csv(file.path(TARGET_dir, "sig_middle_GSEA_results.csv"), row.names = 1)

term = 'regulation of tube size, open tracheal system (GO:0035151)'

target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
norm_exp_1 = plot_heatmap(cds, target_genes)
png(filename = file.path(TARGET_dir, paste0(term, "_dynamic_gene_heatmap.png")), height = 2000, width = 1500, res = 200)
pheatmap(norm_exp_1, cluster_cols = FALSE, cluster_rows = FALSE, main = stringr::str_split(term, ' \\(')[[1]][1], fontsize = 14)
dev.off()

term = 'Golgi vesicle transport (GO:0048193)'
target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
norm_exp_1 = plot_heatmap(cds, target_genes)
png(filename = file.path(TARGET_dir, paste0(term, "_dynamic_gene_heatmap.png")), height = 2000, width = 1000, res = 200)
pheatmap(norm_exp_1, cluster_cols = FALSE, cluster_rows = FALSE, main = stringr::str_split(term, ' \\(')[[1]][1], fontsize = 14)
dev.off()

term = 'chitin-based cuticle development (GO:0040003)'
target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
norm_exp_1 = plot_heatmap(cds, target_genes)
png(filename = file.path(TARGET_dir, paste0(term, "_dynamic_gene_heatmap.png")), height = 2000, width = 1300, res = 200)
pheatmap(norm_exp_1, cluster_cols = FALSE, cluster_rows = FALSE, main = stringr::str_split(term, ' \\(')[[1]][1], fontsize = 14)
dev.off()


##### plot out the transcription factors #####
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "monocle3_no_batch_correct_object.rds"))
early_DE_genes = read.csv("results/v18/early_wt12_enrichment/Trachea/markers_genes.csv", row.names = 1)
late_DE_genes = read.csv("results/v18/wt13_enrichment/Trachea/markers_genes.csv", row.names = 1)
late_DE_genes = late_DE_genes[late_DE_genes$p_val_adj < 0.05 & late_DE_genes$avg_log2FC > 0, ]
early_DE_genes = early_DE_genes[early_DE_genes$p_val_adj < 0.05 & early_DE_genes$avg_log2FC > 0, ]

late_DE_genes$symbol = rownames(late_DE_genes)
early_DE_genes$symbol = rownames(early_DE_genes)

late_DE_genes$type = 'late'
early_DE_genes$type = 'early'
combined_DE_genes = rbind(late_DE_genes, early_DE_genes)
combined_DE_genes[combined_DE_genes$symbol %in% combined_DE_genes$symbol[duplicated(combined_DE_genes$symbol)], 'type'] = 
  paste0(combined_DE_genes[combined_DE_genes$symbol %in% combined_DE_genes$symbol[duplicated(combined_DE_genes$symbol)], 'type'], "_", "both")

TF_tab = read.csv("accessory_data/Drosophila_TFs/all_candidates.csv", sep = '\t')
TF_tab = TF_tab[TF_tab$verdict_DNA_BD != "NO", ]

i_TFs = intersect(TF_tab$symbol, combined_DE_genes$symbol)

combined_DE_genes = combined_DE_genes[combined_DE_genes$symbol %in% i_TFs, ]

p = plot_genes_by_group(cds, markers = i_TFs, norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
  xlab("Cell Types") + 
  coord_flip() + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, 'dynamic_TF.png'), plot = p, width = 30, height = 6)


