library(monocle3)
library(ggplot2)
library(RColorBrewer)
library(dbplyr)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'refined_wt13_early12_trachea')
dir.create(TARGET_dir, recursive = TRUE)

cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "monocle3_no_batch_correct_object.rds"))

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

UMAP_coord$cell_type = cds@colData$cell_type
UMAP_coord[UMAP_coord$cell_type == 'Branching Trachea Cells', 'cell_type'] = 'Tip Cells'

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
  scale_color_brewer(palette = 'Set3')
ggsave(filename = file.path(TARGET_dir, "cluster.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=reorder(batch, pseudotime), y=pseudotime, fill = batch)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set1') + 
  ylab("pseudotime") + 
  xlab("batch")
ggsave(filename = file.path(TARGET_dir, "violin_pseudotime.png"), plot = p, width = 8, height = 6)

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = cell_type)) +
  guides(color=guide_legend(title="")) +
  geom_point() + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set2', breaks=c('Early Trachea Cells', 'Middle Trachea Cells', 'Tip Cells', 'Late Trachea Cells'))
ggsave(filename = file.path(TARGET_dir, "celltypes.png"), plot = p, width = 8, height = 6)

norm_exp = monocle3::normalized_counts(cds)

UMAP_coord$btl = norm_exp['btl', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -btl), y=btl, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("btl normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_btl.png"), plot = p, width = 8, height = 6)

UMAP_coord$Mipp1 = norm_exp['Mipp1', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -Mipp1), y=Mipp1, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("Mipp1 normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_Mipp1.png"), plot = p, width = 8, height = 6)

UMAP_coord$bnl = norm_exp['bnl', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -bnl), y=bnl, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("bnl normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_bnl.png"), plot = p, width = 8, height = 6)

UMAP_coord$pnt = norm_exp['pnt', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -pnt), y=pnt, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("pnt normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_pnt.png"), plot = p, width = 8, height = 6)

UMAP_coord$sty = norm_exp['sty', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -sty), y=sty, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("sty normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_sty.png"), plot = p, width = 8, height = 6)

UMAP_coord$shg = norm_exp['shg', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -shg), y=shg, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("shg normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_shg.png"), plot = p, width = 8, height = 6)

# plot out the tracheal tip cells 
GO_terms = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Branching Trachea Cells_sig_GO_biological_early.csv"), row.names = 1)
GO_terms = GO_terms[order(GO_terms$Adjusted.P.value, decreasing = FALSE), ]
GO_terms = GO_terms[1:40, ]
GO_terms = GO_terms[grepl('axon', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[grepl('wing', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[grepl('glial', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[grepl('eye', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[grepl('nervous', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[grepl('salivary', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[grepl('border follicle', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[grepl('open tracheal system development ', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[grepl('respiratory system development', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[grepl('establishment of blood-brain barrier ', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[grepl('post-embryonic appendage morphogenesis', GO_terms$Term) == FALSE, ]
GO_terms = GO_terms[1:15, ]

GO_terms$log_pval = -log10(GO_terms$Adjusted.P.value)

p = ggplot(data=GO_terms, aes(x=reorder(Term, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[4]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw()
ggsave(filename = file.path(TARGET_dir, "Tip_Cells_EnrichR_results.png"), plot = p, width = 10, height = 4)

# the below are for early cells 
GO_terms = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Early Trachea Cells_sig_GO_biological_early.csv"), row.names = 1)
GO_terms = GO_terms[order(GO_terms$Adjusted.P.value, decreasing = FALSE), ]
GO_terms = GO_terms[1:40, ]
GO_interest = c("mRNA splicing, via spliceosome (GO:0000398)", 
                "chromatin remodeling (GO:0006338)", 
                "regulation of mitotic cell cycle (GO:0007346)", 
                "Notch signaling pathway (GO:0007219)", 
                "tracheal outgrowth, open tracheal system (GO:0007426)")
GO_terms = GO_terms[GO_terms$Term %in% GO_interest, ]
GO_terms$log_pval = -log10(GO_terms$Adjusted.P.value)
p = ggplot(data=GO_terms, aes(x=reorder(Term, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[1]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw()
ggsave(filename = file.path(TARGET_dir, "Early_Cells_EnrichR_results.png"), plot = p, width = 10, height = 4)

# the below are for middle cells 
GO_terms = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Middle Trachea Cells_sig_GO_biological_early.csv"), row.names = 1)
GO_terms = GO_terms[order(GO_terms$Adjusted.P.value, decreasing = FALSE), ]
GO_terms = GO_terms[1:40, ]
GO_interest = c("regulation of tube size, open tracheal system (GO:0035151)", 
                "cell-cell junction assembly (GO:0007043)", 
                "Golgi vesicle transport (GO:0048193)", 
                "dorsal closure (GO:0007391)", 
                "Ras protein signal transduction (GO:0007265)")
GO_terms = GO_terms[GO_terms$Term %in% GO_interest, ]
GO_terms$log_pval = -log10(GO_terms$Adjusted.P.value)
p = ggplot(data=GO_terms, aes(x=reorder(Term, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[3]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw()
ggsave(filename = file.path(TARGET_dir, "Middle_Cells_EnrichR_results.png"), plot = p, width = 10, height = 4)

# the below are for late cells 
GO_terms = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Late Trachea Cells_sig_GO_biological_early.csv"), row.names = 1)
GO_terms = GO_terms[order(GO_terms$Adjusted.P.value, decreasing = FALSE), ]
GO_terms = GO_terms[1:40, ]
GO_interest = c("cytoplasmic translation (GO:0002181)", 
                "ribosome assembly (GO:0042255)", 
                "cellular protein metabolic process (GO:0044267)")
GO_terms = GO_terms[GO_terms$Term %in% GO_interest, ]
GO_terms$log_pval = -log10(GO_terms$Adjusted.P.value)
p = ggplot(data=GO_terms, aes(x=reorder(Term, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 4, 'Set2')[2]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw()
ggsave(filename = file.path(TARGET_dir, "Late_Cells_EnrichR_results.png"), plot = p, width = 10, height = 4)

#############################
# plot out the Golgi Vesicle gene expression 
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "monocle3_no_batch_correct_object.rds"))
cds = cds[, monocle3::clusters(cds) != 9] # remove the branch from the main trajectory
term = 'Golgi vesicle transport (GO:0048193)'
GO_terms = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Middle Trachea Cells_sig_GO_biological_early.csv"), row.names = 1)

target_genes = GO_terms[GO_terms$Term == term, 'Genes']
target_genes = stringr::str_split(target_genes, pattern = ";")[[1]]

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
  geom_line(color = RColorBrewer::brewer.pal(n = 4, 'Set2')[3]) + theme_bw() 
ggsave(file.path(TARGET_dir, paste0(term, "_dynamic_gene_line_avg.png")), plot = p, width = 8, height = 5)

# plot out the Cytoplasmic translation
term = 'cytoplasmic translation (GO:0002181)'
GO_terms = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Late Trachea Cells_sig_GO_biological_early.csv"), row.names = 1)
target_genes = GO_terms[GO_terms$Term == term, 'Genes']
target_genes = stringr::str_split(target_genes, pattern = ";")[[1]]

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


##############################################
# this is plot out the dynamically expressed TFs 
ATres = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", 'raw_associationTest.rds'))
ATres = ATres[!is.na(ATres$pvalue), ]
ATres$adj_p = p.adjust(ATres$pvalue, method = 'fdr')
ATres = ATres[ATres$adj_p < 0.05, ]
pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')

TF_tab = read.csv("accessory_data/Drosophila_TFs/all_candidates.csv", sep = '\t')
TF_tab = TF_tab[TF_tab$verdict_DNA_BD != "NO", ]
i_TFs = intersect(TF_tab$symbol, rownames(ATres))

target_genes = i_TFs
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
png(filename = file.path(TARGET_dir, paste0("TF_dynamic_gene_heatmap.png")), height = 2000, width = 1000, res = 200)
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
rank_sum_test = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", 'rank_sum_test.csv'), row.names = 1)
rank_sum_test = rank_sum_test[rank_sum_test$padj < 0.05, ]
rank_sum_test = rank_sum_test[rank_sum_test$logFC > 0.1, ]
unique(plot_df$gene[plot_df$gene %in% rank_sum_test$feature])
#[1] "aop"        "cg"         "psq"        "grh"        "Ssb-c31a"   "Su(var)205" "HmgZ"       "Dsp1"      
#[9] "cic"        "rib"        "Pep"        "Eip75B"     "noc"        "bowl"       "grn"        "ct"        
#[17] "CrebA"      "Xbp1"       "ttk"        "lola"       "EcR"        "Hr4"        "luna"     

#sg_TFs = intersect(pathway_list[["salivary gland development (GO:0007431)"]], TF_tab$symbol)
scaled_exp = scaled_exp[rownames(scaled_exp) %in% unique(plot_df$gene[plot_df$gene %in% rank_sum_test$feature]), ]
plot_df = convert_line_plot(scaled_exp)

plot_df$pseudotime = (plot_df$pseudotime - min(plot_df$pseudotime)) / max(plot_df$pseudotime)

# interesting genes 
# genes in development of trachea 
Trachea_genes = pathway_list[["open tracheal system development (GO:0007424)"]]

sub_plot_df = plot_df
sub_plot_df = sub_plot_df[sub_plot_df$gene %in% Trachea_genes, ]
p<-ggplot(sub_plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
  xlab("pseudotime") + 
  ylab("average expression") +
  geom_line(aes(color = gene)) + theme_bw() 
ggsave(file.path(TARGET_dir, "list_A_dynamic_gene_line_avg.png"), plot = p, width = 8, height = 5)














