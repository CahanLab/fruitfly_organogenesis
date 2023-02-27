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
p = ggplot(gsea_results, aes(x = reorder(pathway, -padj), y = NES)) + 
  geom_point(aes(size = GeneRatio, color = padj)) +
  scale_size_continuous(range = c(4,8)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab('Normalized Enrichment Scores') +
  xlab("GO terms") +
  ggtitle("Stage 13-16: Salivary Gland Enrichment") + 
  coord_flip()
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
p = ggplot(gsea_results, aes(x = reorder(pathway, -padj), y = NES)) + 
  geom_point(aes(size = GeneRatio, color = padj)) +
  scale_size_continuous(range = c(4,8)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab('Normalized Enrichment Scores') +
  xlab("GO terms") +
  ggtitle("Stage 13-16: Trachea Enrichment") + 
  coord_flip()
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
p = ggplot(gsea_results, aes(x = reorder(pathway, -padj), y = NES)) + 
  geom_point(aes(size = GeneRatio, color = padj)) +
  scale_size_continuous(range = c(4,8)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab('Normalized Enrichment Scores') +
  xlab("GO terms") +
  ggtitle("Stage 13-16: Germ Cells Enrichment") + 
  coord_flip()
ggsave(filename = file.path(TARGET_dir, "stage13-16_gc.png"), width = 12, height = 10)
