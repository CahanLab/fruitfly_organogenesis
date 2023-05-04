library(ggplot2)
library(RColorBrewer)
library(ggdendroplot)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'wt13_GSEA')

gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Salivary Gland/gsea_results_wt.csv"), row.names = 1)
gsea_results = gsea_results[gsea_results$padj < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]
gsea_results$num_leadingEdges = NULL
for(temp_index in rownames(gsea_results)) {
  leading_edge = eval(parse(text = gsea_results[temp_index, 'leadingEdge']))
  gsea_results[temp_index, 'num_leadingEdges'] = length(leading_edge)
}
gsea_results$GeneRatio = gsea_results$num_leadingEdges / gsea_results$size
gsea_results$NES = as.numeric(gsea_results$NES)
gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
gsea_results$pathway = stringr::str_replace_all(gsea_results$pathway, "\\(", "\n\\(\\")

p = ggplot(gsea_results, aes(x = reorder(pathway, -padj), y = NES)) + 
  geom_point(aes(size = GeneRatio, color = padj)) +
  scale_size_continuous(range = c(4,8)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab('Normalized Enrichment Scores') +
  xlab("GO terms") +
  ggtitle("Stage 13-16: Salivary Gland Enrichment") + 
  coord_flip() + 
  theme(text = element_text(size = 24), plot.title.position = "plot")
ggsave(filename = file.path(TARGET_dir, "stage13-16_sg.png"), width = 12, height = 10)

# this is to plot out trachea 
gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Trachea/gsea_results_wt.csv"), row.names = 1)
gsea_results = gsea_results[gsea_results$padj < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]
gsea_results$num_leadingEdges = NULL
for(temp_index in rownames(gsea_results)) {
  leading_edge = eval(parse(text = gsea_results[temp_index, 'leadingEdge']))
  gsea_results[temp_index, 'num_leadingEdges'] = length(leading_edge)
}
gsea_results$GeneRatio = gsea_results$num_leadingEdges / gsea_results$size
gsea_results$NES = as.numeric(gsea_results$NES)
gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
gsea_results = gsea_results[1:10, ]
gsea_results$pathway = stringr::str_replace_all(gsea_results$pathway, "\\(", "\n\\(\\")

p = ggplot(gsea_results, aes(x = reorder(pathway, -padj), y = NES)) + 
  geom_point(aes(size = GeneRatio, color = padj)) +
  scale_size_continuous(range = c(4,8)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab('Normalized Enrichment Scores') +
  xlab("GO terms") +
  ggtitle("Stage 13-16: Trachea Enrichment") + 
  coord_flip() + 
  theme(text = element_text(size = 24), plot.title.position = "plot")

ggsave(filename = file.path(TARGET_dir, "stage13-16_trachea.png"), width = 12, height = 10)

# this is to plot out germ cells 
gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Germ Cells/gsea_results_wt.csv"), row.names = 1)
gsea_results = gsea_results[gsea_results$padj < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]
gsea_results$num_leadingEdges = NULL
for(temp_index in rownames(gsea_results)) {
  leading_edge = eval(parse(text = gsea_results[temp_index, 'leadingEdge']))
  gsea_results[temp_index, 'num_leadingEdges'] = length(leading_edge)
}
gsea_results$GeneRatio = gsea_results$num_leadingEdges / gsea_results$size
gsea_results$NES = as.numeric(gsea_results$NES)
gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
gsea_results = gsea_results[1:10, ]
gsea_results[gsea_results$pathway == 'proteasome-mediated ubiquitin-dependent protein catabolic process (GO:0043161)', 'pathway'] = 'proteasome-mediated ubiquitin-dependent \n protein catabolic process (GO:0043161)'
gsea_results$pathway = stringr::str_replace_all(gsea_results$pathway, "\\(", "\n\\(\\")

p = ggplot(gsea_results, aes(x = reorder(pathway, -padj), y = NES)) + 
  geom_point(aes(size = GeneRatio, color = padj)) +
  scale_size_continuous(range = c(4,8)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab('Normalized Enrichment Scores') +
  xlab("GO terms") +
  ggtitle("Stage 13-16: Germ Cells Enrichment") + 
  coord_flip() + 
  theme(text = element_text(size = 24), plot.title.position = "plot")

ggsave(filename = file.path(TARGET_dir, "stage13-16_gc.png"), width = 16, height = 10)


##### this is to plot out a individual barplot #####
gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Salivary Gland/gsea_results_wt.csv"), row.names = 1)
gsea_results = gsea_results[gsea_results$padj < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]
gsea_results$log_padj = -log10(gsea_results$padj)

interesting_cat = c("Golgi vesicle transport (GO:0048193)", 
                    "cytoplasmic translation (GO:0002181)", 
                    "exocrine system development (GO:0035272)")
sub_gsea_results = gsea_results[gsea_results$pathway %in% interesting_cat, ]
sub_gsea_results$pathway = stringr::str_replace_all(sub_gsea_results$pathway, "\\(", "\n\\(\\")

big_plot_df = data.frame("pathway" = sub_gsea_results$pathway, 
                         'log10_adj' = sub_gsea_results$log_padj, 
                         'celltype' = 'salivary \n gland')

gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Trachea/gsea_results_wt.csv"), row.names = 1)
gsea_results = gsea_results[gsea_results$padj < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]
gsea_results$log_padj = -log10(gsea_results$padj)

interesting_cat = c("open tracheal system development (GO:0007424)", 
                    "regulation of tube size, open tracheal system (GO:0035151)", 
                    "cell-cell junction assembly (GO:0007043)", 
                    "Golgi vesicle transport (GO:0048193)", 
                    "liquid clearance, open tracheal system (GO:0035002)", 
                    "epithelial cell migration, open tracheal system (GO:0007427)", 
                    "establishment of epithelial cell apical/basal polarity (GO:0045198)")
sub_gsea_results = gsea_results[gsea_results$pathway %in% interesting_cat, ]
sub_gsea_results$pathway = stringr::str_replace_all(sub_gsea_results$pathway, "\\(", "\n\\(\\")

temp_plot_df = data.frame("pathway" = sub_gsea_results$pathway, 
                          'log10_adj' = sub_gsea_results$log_padj, 
                          'celltype' = 'trachea')
big_plot_df = rbind(big_plot_df, temp_plot_df)

gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Germ Cells/gsea_results_wt.csv"), row.names = 1)
gsea_results = gsea_results[gsea_results$padj < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]
gsea_results$log_padj = -log10(gsea_results$padj)

interesting_cat = c("proteasomal protein catabolic process (GO:0010498)", 
                    "ubiquitin-dependent protein catabolic process (GO:0006511)")
sub_gsea_results = gsea_results[gsea_results$pathway %in% interesting_cat, ]
sub_gsea_results$pathway = stringr::str_replace_all(sub_gsea_results$pathway, "\\(", "\n\\(\\")

temp_plot_df = data.frame("pathway" = sub_gsea_results$pathway, 
                          'log10_adj' = sub_gsea_results$log_padj, 
                          'celltype' = 'germ \n cells')
big_plot_df = rbind(big_plot_df, temp_plot_df)
big_plot_df$pathway = stringr::str_remove_all(big_plot_df$pathway, " \\\n.*")
#big_plot_df[big_plot_df$pathway == "cell morphogenesis involved in differentiation", "pathway"] = "cell morphogenesis \n involved in differentiation"
#big_plot_df[big_plot_df$pathway == "establishment or maintenance of apical/basal cell polarity", "pathway"] = "establishment or maintenance of \n apical/basal cell polarity"
#big_plot_df[big_plot_df$pathway == "proteasomal ubiquitin-independent protein catabolic process", "pathway"] = "proteasomal ubiquitin-independent \n protein catabolic process"

big_plot_df$celltype = factor(big_plot_df$celltype, levels = c("salivary \n gland", "trachea", "germ \n cells"))
p <- ggplot(data = big_plot_df, aes(y = reorder(pathway, log10_adj), x = log10_adj, fill = celltype)) +
  geom_bar(stat="identity") +
  #scale.func(range = c(0, 100), limits = c(scale.min, scale.max)) +
  #guides(size = guide_legend(title = 'Percent Expressed')) +
  #guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
  #scale_colour_viridis_c() + 
  labs(
    x = '-log10 adjusted p-value',
    y = ''
  ) + 
  theme_classic()  + 
  facet_grid(
    rows = vars(celltype),
    scales = "free_y",
    space = "free_y",
    switch = "x"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank(), 
    text = element_text(size = 25), 
    legend.position="none"
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=0, vjust = 1, hjust=1)) +
  ggtitle("Stage 13-16 Embryos")  
ggsave(filename = file.path(TARGET_dir, "GO_terms_embryo.png"), width = 12, height = 8)


