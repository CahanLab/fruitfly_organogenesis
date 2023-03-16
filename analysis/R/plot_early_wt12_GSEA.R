library(ggplot2)
library(RColorBrewer)
library(ggdendroplot)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'early_wt12_GSEA')

gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Salivary Gland/gsea_results_wt.csv"), row.names = 1)
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
  ggtitle("Stage 10-12: Salivary Gland Enrichment") + 
  coord_flip() + 
  theme(text = element_text(size = 24), plot.title.position = "plot")

ggsave(filename = file.path(TARGET_dir, "stage10-12_sg.png"), width = 12, height = 10)

# this is to plot out trachea 
gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Trachea/gsea_results_wt.csv"), row.names = 1)
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
  ggtitle("Stage 10-12: Trachea Enrichment") + 
  coord_flip() + 
  theme(text = element_text(size = 24), plot.title.position = "plot")

ggsave(filename = file.path(TARGET_dir, "stage10-12_trachea.png"), width = 12, height = 10)

# this is to plot out germ cells 
gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Germ Cells/gsea_results_wt.csv"), row.names = 1)
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
  ggtitle("Stage 10-12: Germ Cells Enrichment") + 
  coord_flip() + 
  theme(text = element_text(size = 24), plot.title.position = "plot")

ggsave(filename = file.path(TARGET_dir, "stage10-12_gc.png"), width = 12, height = 10)
