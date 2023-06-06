library(monocle3)
library(ggplot2)
library(RColorBrewer)
library(dbplyr)
library(ggsignif)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'refined_wt13_early12_germ')
dir.create(TARGET_dir, recursive = TRUE)

rank_sum = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ", "rank_sum_test.csv"), row.names = 1)
rank_sum = rank_sum[rank_sum$logFC > 0, ]
rank_sum[rank_sum$group == '1', 'group'] = 'Unknown 1'
rank_sum[rank_sum$group == '3', 'group'] = 'Unknown 2'
rank_sum[rank_sum$group == '5', 'group'] = 'Early Germ Cells'
rank_sum[rank_sum$group == '6', 'group'] = 'Interm. Germ Cells 1'
rank_sum[rank_sum$group == '2', 'group'] = 'Interm. Germ Cells 2'
rank_sum[rank_sum$group == '4', 'group'] = 'Late Germ Cells'

write.csv(rank_sum, file = file.path(TARGET_dir, 'DE_genes.csv'))

cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ", "monocle3_no_batch_correct_object.rds"))
plot_cells(cds, genes = 'BigH1', cell_size = 1) 
plot_cells(cds, genes = c("lncRNA:roX1", "lncRNA:roX2"), cell_size = 1, scale_to_range = FALSE, show_trajectory_graph = FALSE) 
plot_cells(cds, genes = c("CG6701", "Pp2C1", 'ovo', 'otu', 'Sxl'), cell_size = 1, scale_to_range = FALSE, show_trajectory_graph = FALSE) 
plot_cells(cds, genes = c("FDY"), cell_size = 1, scale_to_range = FALSE, show_trajectory_graph = FALSE) 

# this is to plot out the dead cells 
plot_cells(cds, genes = c('nos', 'wun2', 'Lsd-1', 'Lsd-2'), cell_size = 1, show_trajectory_graph = FALSE) 

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

UMAP_coord$cell_type = UMAP_coord$clusters
UMAP_coord[UMAP_coord$clusters == '1', 'cell_type'] = 'Unknown 1'
UMAP_coord[UMAP_coord$clusters == '3', 'cell_type'] = 'Unknown 2'
UMAP_coord[UMAP_coord$clusters == '5', 'cell_type'] = 'Early Germ Cells'
UMAP_coord[UMAP_coord$clusters == '6', 'cell_type'] = 'Interm. Germ Cells 1'
UMAP_coord[UMAP_coord$clusters == '2', 'cell_type'] = 'Interm. Germ Cells 2'
UMAP_coord[UMAP_coord$clusters == '4', 'cell_type'] = 'Late Germ Cells'


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

UMAP_coord$cell_type <- factor(UMAP_coord$cell_type, levels = c("Early Germ Cells", 
                                                                "Interm. Germ Cells 1", 
                                                                "Interm. Germ Cells 2", 
                                                                "Late Germ Cells", 
                                                                "Unknown 1", 
                                                                "Unknown 2"))

p = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = cell_type)) +
  guides(color=guide_legend(title="")) +
  geom_point() + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set2') + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "celltypes.png"), plot = p, width = 8, height = 6)

##### plotting the gene sets results #####
norm_exp = monocle3::normalized_counts(cds)
gene_interest_list = c('eya', 'wg', 'tkv', 'Dl', 'shg', 'Prosalpha5', 'Sod1', 'BigH1', 'bru1', 'dsx', 'dhd', 'vas', 'pum', 
                       'nos', 'bam', 'osk', 'pgc', 'gcl', 'nos', 'dpp', 'N', 'Pomp')
for(gene_interest in gene_interest_list) {
  UMAP_coord$gene_exp = norm_exp[gene_interest, ]
  p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -gene_exp), y=gene_exp, fill = cell_type)) + 
    geom_violin() +
    guides(fill=guide_legend(title="")) +
    geom_boxplot(width=0.1) +
    theme_minimal() +
    scale_fill_brewer(palette = 'Set2') + 
    ylab(paste0(gene_interest, " normalized expression")) + 
    xlab("cell type") + 
    theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(filename = file.path(TARGET_dir, paste0("violin_", gene_interest, ".png")), plot = p, width = 8, height = 6)
}

gene_interest_list = c('bam', 'bgcn', 'zpg', 'stg', 'esg', 'kmg', 'Rbp4', 'Pp2C1', 'CG6701', 'Parp16', 'Glut3')
for(gene_interest in gene_interest_list) {
  UMAP_coord$gene_exp = norm_exp[gene_interest, ]
  p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -gene_exp), y=gene_exp, fill = cell_type)) + 
    geom_violin() +
    guides(fill=guide_legend(title="")) +
    geom_boxplot(width=0.1) +
    theme_minimal() +
    scale_fill_brewer(palette = 'Set2') + 
    ylab(paste0(gene_interest, " normalized expression")) + 
    xlab("cell type") + 
    theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(filename = file.path(TARGET_dir, paste0("violin_", gene_interest, ".png")), plot = p, width = 8, height = 6)
}

gene_interest_list = c("lncRNA:roX1", "lncRNA:roX2")
for(gene_interest in gene_interest_list) {
  UMAP_coord$gene_exp = norm_exp[gene_interest, ]
  p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -gene_exp), y=gene_exp, fill = cell_type)) + 
    geom_violin() +
    guides(fill=guide_legend(title="")) +
    geom_boxplot(width=0.1) +
    theme_minimal() +
    scale_fill_brewer(palette = 'Set2') + 
    ylab(paste0(gene_interest, " normalized expression")) + 
    xlab("cell type") + 
    theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(filename = file.path(TARGET_dir, paste0("violin_", gene_interest, ".png")), plot = p, width = 8, height = 6)
}
###### plot out the GSEA results #####
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "4_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_late_GSEA_results.csv"))

focus_gsea = c("proteasomal protein catabolic process (GO:0010498)", 
               "hydrogen ion transmembrane transport (GO:1902600)", 
               'mitochondrial transport (GO:0006839)', 
               'ATP biosynthetic process (GO:0006754)', 
               'ATP hydrolysis coupled cation transmembrane transport (GO:0099132)')

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results[sub_GSEA_results$pathway == "ATP hydrolysis coupled cation transmembrane transport (GO:0099132)", 'pathway'] = 'ATP hydrolysis coupled \n cation transmembrane transport (GO:0099132)'
sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 6, 'Set2')[4]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "late_Cells_GSEA_results.png"), plot = p, width = 10, height = 4.7)

# middle stage 2
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "2_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_middle_2_GSEA_results.csv"))

focus_gsea = c("proteasomal protein catabolic process (GO:0010498)", 
               "female gamete generation (GO:0007292)", 
               "negative regulation of translation (GO:0017148)", 
               'pole cell migration (GO:0007280)', 
               'male gamete generation (GO:0048232)', 
               'ATP hydrolysis coupled cation transmembrane transport (GO:0099132)')

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results[sub_GSEA_results$pathway == "ATP hydrolysis coupled cation transmembrane transport (GO:0099132)", 'pathway'] = 'ATP hydrolysis coupled \n cation transmembrane transport (GO:0099132)'

sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 6, 'Set2')[3]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw() +   
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "middle_2_Cells_GSEA_results.png"), plot = p, width = 10, height = 5.5)

term = "negative regulation of translation (GO:0017148)"
target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
genes_translation = target_genes

term = "negative regulation of transcription, DNA-templated (GO:0045892)"
target_genes = GSEA_results[GSEA_results$pathway == term, 'leadingEdge']
target_genes = eval(parse(text = target_genes))
genes_transcription = target_genes

length(intersect(genes_transcription, genes_translation))

# middle stage 1
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "6_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_middle_1_GSEA_results.csv"))

focus_gsea = c("mRNA splicing, via spliceosome (GO:0000398)",
               "cytoskeleton-dependent cytokinesis (GO:0061640)", 
               "mitotic cytokinesis (GO:0000281)", 
               "cytoplasmic translation (GO:0002181)", 
               "protein catabolic process (GO:0030163)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 6, 'Set2')[2]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "middle_1_Cells_GSEA_results.png"), plot = p, width = 9.5, height = 4)

# early stage
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "5_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_early_GSEA_results.csv"))

focus_gsea = c("mRNA splicing, via spliceosome (GO:0000398)",
               "regulation of gene expression (GO:0010468)", 
               "Notch signaling pathway (GO:0007219)", 
               "regulation of macromolecule metabolic process (GO:0060255)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 6, 'Set2')[1]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "early_Cells_GSEA_results.png"), plot = p, width = 10, height = 4)

# unknown 1 stage
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "1_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_unknown_1_GSEA_results.csv"))

focus_gsea = c("cytoplasmic translation (GO:0002181)",
               "ribosome assembly (GO:0042255)", 
               "ribosomal small subunit biogenesis (GO:0042274)", 
               "ribosomal large subunit biogenesis (GO:0042273)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 8, 'Set2')[5]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "unknown_1_Cells_GSEA_results.png"), plot = p, width = 10, height = 4)

# unknown 2 stage
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "3_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_unknown_2_GSEA_results.csv"))

focus_gsea = c("positive regulation of transcription, DNA-templated (GO:0045893)",
               "Wnt signaling pathway (GO:0016055)", 
               "cell-cell adhesion via plasma-membrane adhesion molecules (GO:0098742)",
               "imaginal disc-derived wing morphogenesis (GO:0007476)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)
sub_GSEA_results[sub_GSEA_results$pathway == "cell-cell adhesion via plasma-membrane adhesion molecules (GO:0098742)", 'pathway'] = "cell-cell adhesion via plasma-membrane \n adhesion molecules (GO:0098742)"
sub_GSEA_results$pathway = stringr::str_replace_all(sub_GSEA_results$pathway, "\\(", "\n\\(\\")

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 8, 'Set2')[6]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, "unknown_2_Cells_GSEA_results.png"), plot = p, width = 10, height = 4)

##### make some dotplots #####
cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ", "monocle3_no_batch_correct_object.rds"))
cds@colData$cell_type = NA
cds@colData[clusters(cds) == '1', 'cell_type'] = 'Unknown 1'
cds@colData[clusters(cds) == '3', 'cell_type'] = 'Unknown 2'
cds@colData[clusters(cds) == '5', 'cell_type'] = 'Early Germ Cells'
cds@colData[clusters(cds) == '6', 'cell_type'] = 'Interm. Germ Cells 1'
cds@colData[clusters(cds) == '2', 'cell_type'] = 'Interm. Germ Cells 2'
cds@colData[clusters(cds) == '4', 'cell_type'] = 'Late Germ Cells'

p = plot_genes_by_group(cds, markers = c('nos', 'gcl', 'pgc', 'vas'), norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
  xlab("Cell Types") + 
  coord_flip() + 
  scale_x_discrete(limits = c('Unknown 2', 'Unknown 1', 'Late Germ Cells', 'Interm. Germ Cells 2', 'Interm. Germ Cells 1', 'Early Germ Cells')) + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, 'germ_stem_cell_markers.png'), plot = p, width = 10, height = 5)

# geensets for germ-line stem cell asymmetric division
set_genes = c("Ote", "eff", "aurB", "pelo", "twin", "eIF4A")
p = plot_genes_by_group(cds, markers = set_genes, norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
  xlab("Cell Types") + 
  coord_flip() + 
  scale_x_discrete(limits = c('Unknown 2', 'Unknown 1', 'Late Germ Cells', 'Interm. Germ Cells 2', 'Interm. Germ Cells 1', 'Early Germ Cells')) + 
  theme(text = element_text(size = 24))

# somatic female genes
set_genes = c("Sxl", "tra")
p = plot_genes_by_group(cds, markers = set_genes, norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
  xlab("Cell Types") + 
  coord_flip() + 
  scale_x_discrete(limits = c('Unknown 2', 'Unknown 1', 'Late Germ Cells', 'Interm. Germ Cells 2', 'Interm. Germ Cells 1', 'Early Germ Cells')) + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, 'female_somatic_markers.png'), plot = p, width = 8, height = 5)

# male germ genes 
set_genes = c("dpa", "Rs1", "Mcm5", "CG9253", "CG6693", "nclb", "Klp61F", "CG6701", "Pp2C1")
p = plot_genes_by_group(cds, markers = set_genes, norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
  xlab("Cell Types") + 
  coord_flip() + 
  scale_x_discrete(limits = c('Unknown 2', 'Unknown 1', 'Late Germ Cells', 'Interm. Germ Cells 2', 'Interm. Germ Cells 1', 'Early Germ Cells')) + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, 'male_germ_markers.png'), plot = p, width = 12, height = 5)

# female germ genes 
set_genes = c("ovo", "otu", "Sxl")
p = plot_genes_by_group(cds, markers = set_genes, norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
  xlab("Cell Types") + 
  coord_flip() + 
  scale_x_discrete(limits = c('Unknown 2', 'Unknown 1', 'Late Germ Cells', 'Interm. Germ Cells 2', 'Interm. Germ Cells 1', 'Early Germ Cells')) + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, 'female_germ_markers.png'), plot = p, width = 8, height = 5.5)

# death genes
set_genes = c("wun2", "Lsd-1", 'Lsd-2')
p = plot_genes_by_group(cds, markers = set_genes, norm_method = 'log', group_cells_by = 'cell_type', ordering_type = 'none') + 
  xlab("Cell Types") + 
  coord_flip() + 
  scale_x_discrete(limits = c('Unknown 2', 'Unknown 1', 'Late Germ Cells', 'Middle Germ Cells 2', 'Middle Germ Cells 1', 'Early Germ Cells')) + 
  theme(text = element_text(size = 24))
ggsave(filename = file.path(TARGET_dir, 'female_germ_markers.png'), plot = p, width = 8, height = 5)

##### test the sequencing depth #####
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

UMAP_coord$cell_type = UMAP_coord$clusters
UMAP_coord[UMAP_coord$clusters == '1', 'cell_type'] = 'Unknown 1'
UMAP_coord[UMAP_coord$clusters == '3', 'cell_type'] = 'Unknown 2'
UMAP_coord[UMAP_coord$clusters == '5', 'cell_type'] = 'Early Germ Cells'
UMAP_coord[UMAP_coord$clusters == '6', 'cell_type'] = 'Interm. Germ Cells 1'
UMAP_coord[UMAP_coord$clusters == '2', 'cell_type'] = 'Interm. Germ Cells 2'
UMAP_coord[UMAP_coord$clusters == '4', 'cell_type'] = 'Late Germ Cells'

UMAP_coord$clusters <- factor(UMAP_coord$clusters, levels = c("Early Germ Cells", 
                                                                "Interm. Germ Cells 1", 
                                                                "Interm. Germ Cells 2", 
                                                                "Late Germ Cells", 
                                                                "Unknown 1", 
                                                                "Unknown 2"))
UMAP_coord$log10_ngenes = log(cds@colData$nCount_RNA)
UMAP_coord$ngenes = cds@colData$nCount_RNA

p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -log10_ngenes), y=log10_ngenes, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("log(UMI Counts)") + 
  xlab("cell type") + 
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_signif(comparisons = list(c("Unknown 2", "Late Germ Cells"), 
                                 c("Unknown 2", "Interm. Germ Cells 1"), 
                                 c("Unknown 2", "Interm. Germ Cells 2"), 
                                 c('Late Germ Cells', 'Unknown 1'), 
                                 c('Interm. Germ Cells 1', 'Unknown 1'),
                                 c('Interm. Germ Cells 2', 'Unknown 1'),
                                 c('Unknown 1', 'Early Germ Cells'), 
                                 c('Unknown 2', 'Early Germ Cells')), map_signif_level=TRUE, step_increase = 0.075)
ggsave(filename = file.path(TARGET_dir, 'read_depth.png'), plot = p, width = 10, height = 5)

median(UMAP_coord[UMAP_coord$cell_type == 'Unknown 1', 'ngenes'])
median(UMAP_coord[UMAP_coord$cell_type == 'Unknown 2', 'ngenes'])

##### check ribosomal biogenesis genes ##### 
ribo_genesis_df = read.csv("accessory_data/Sanchez_et_al_2015/differentiation_enrichment.csv")
# load in the necessary code to cell type 
gtf_last_field = unique(read.table("../quantification/reference_genome_info/dmel-all-r6.33.gtf.gz", sep = "\t")[["V9"]])
gene_metadata = data.frame(
  flybase_id =
    gtf_last_field %>%
    strsplit(";") %>%
    sapply(extract2, 1) %>%
    gsub("^.*gene_id ", "", .),
  symbol =
    gtf_last_field %>%
    strsplit(";") %>%
    sapply(extract2, 2) %>%
    gsub("^.*gene_symbol ", "", .)
)
gene_metadata = gene_metadata %>% dplyr::distinct()
gene_converter = as.vector(gene_metadata$symbol)
names(gene_converter) = as.vector(gene_metadata$flybase_id)

ribo_genes = ribo_genesis_df[ribo_genesis_df$GO.term == 'translation', "Genes..Flybase.ID."]
ribo_genes = stringr::str_split(ribo_genes, ", ")[[1]]
ribo_genes = stringr::str_replace_all(ribo_genes, "GN", "gn")
gene_converter[ribo_genes]
