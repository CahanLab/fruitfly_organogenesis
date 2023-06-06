# make plots for GSEA results for early salivary gland relative to rest of the embryonic cells 
# Fig 2D
TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'early_wt12_GSEA')
dir.create(TARGET_dir, recursive = TRUE)

##### this is to plot out a individual barplot for selected GSEA #####
gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Salivary Gland/gsea_results_wt.csv"), row.names = 1)
gsea_results = gsea_results[gsea_results$padj < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]
gsea_results$log_padj = -log10(gsea_results$padj)

interesting_cat = c("Golgi vesicle transport (GO:0048193)", 
                    "exocrine system development (GO:0035272)", 
                    "salivary gland morphogenesis (GO:0007435)", 
                    "apical junction assembly (GO:0043297)")
sub_gsea_results = gsea_results[gsea_results$pathway %in% interesting_cat, ]
sub_gsea_results$pathway = stringr::str_replace_all(sub_gsea_results$pathway, "\\(", "\n\\(\\")

big_plot_df = data.frame("pathway" = sub_gsea_results$pathway, 
                         'log10_adj' = sub_gsea_results$log_padj, 
                         'celltype' = 'salivary \n gland')

gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Trachea/gsea_results_wt.csv"), row.names = 1)
gsea_results = gsea_results[gsea_results$padj < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]
gsea_results$log_padj = -log10(gsea_results$padj)

interesting_cat = c("open tracheal system development (GO:0007424)", 
                    "regulation of tube size, open tracheal system (GO:0035151)", 
                    "cell morphogenesis involved in differentiation (GO:0000904)", 
                    "apical junction assembly (GO:0043297)", 
                    "regulation of developmental growth (GO:0048638)", 
                    "establishment or maintenance of apical/basal cell polarity (GO:0035088)")
sub_gsea_results = gsea_results[gsea_results$pathway %in% interesting_cat, ]
sub_gsea_results$pathway = stringr::str_replace_all(sub_gsea_results$pathway, "\\(", "\n\\(\\")

temp_plot_df = data.frame("pathway" = sub_gsea_results$pathway, 
                          'log10_adj' = sub_gsea_results$log_padj, 
                          'celltype' = 'trachea')
big_plot_df = rbind(big_plot_df, temp_plot_df)

gsea_results = read.csv(file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Germ Cells/gsea_results_wt.csv"), row.names = 1)
gsea_results = gsea_results[gsea_results$padj < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]
gsea_results$log_padj = -log10(gsea_results$padj)

interesting_cat = c("proteasomal ubiquitin-independent protein catabolic process (GO:0010499)", 
                    "ubiquitin-dependent protein catabolic process (GO:0006511)", 
                    "ATP hydrolysis coupled cation transmembrane transport (GO:0099132)")
sub_gsea_results = gsea_results[gsea_results$pathway %in% interesting_cat, ]
sub_gsea_results$pathway = stringr::str_replace_all(sub_gsea_results$pathway, "\\(", "\n\\(\\")

temp_plot_df = data.frame("pathway" = sub_gsea_results$pathway, 
                          'log10_adj' = sub_gsea_results$log_padj, 
                          'celltype' = 'germ \n cells')
big_plot_df = rbind(big_plot_df, temp_plot_df)
big_plot_df$pathway = stringr::str_remove_all(big_plot_df$pathway, " \\\n.*")
big_plot_df[big_plot_df$pathway == "establishment or maintenance of apical/basal cell polarity", "pathway"] = "establishment or maintenance of \n apical/basal cell polarity"
big_plot_df[big_plot_df$pathway == "proteasomal ubiquitin-independent protein catabolic process", "pathway"] = "proteasomal ubiquitin-independent \n protein catabolic process"

big_plot_df$celltype = factor(big_plot_df$celltype, levels = c("salivary \n gland", "trachea", "germ \n cells"))
p <- ggplot(data = big_plot_df, aes(y = reorder(pathway, log10_adj), x = log10_adj, fill = celltype)) +
  geom_bar(stat="identity") +
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
  ggtitle("Stage 10-12 Embryos")  
ggsave(filename = file.path(TARGET_dir, "GO_terms_embryo.png"), width = 12, height = 8)

