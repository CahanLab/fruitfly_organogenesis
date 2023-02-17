library(monocle3)
library(ggplot2)
library(RColorBrewer)
library(dbplyr)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'refined_wt13_early12_germ')
dir.create(TARGET_dir, recursive = TRUE)

rank_sum = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ", "rank_sum_test.csv"), row.names = 1)
rank_sum = rank_sum[rank_sum$logFC > 0, ]
rank_sum[rank_sum$group == '1', 'group'] = 'Unknown 1'
rank_sum[rank_sum$group == '3', 'group'] = 'Unknown 2'
rank_sum[rank_sum$group == '5', 'group'] = 'Early Germ Cells'
rank_sum[rank_sum$group == '6', 'group'] = 'Middle Germ Cells 1'
rank_sum[rank_sum$group == '2', 'group'] = 'Middle Germ Cells 2'
rank_sum[rank_sum$group == '4', 'group'] = 'Late Germ Cells 1'

write.csv(rank_sum, file = file.path(TARGET_dir, 'DE_genes.csv'))

cds = readRDS(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ", "monocle3_no_batch_correct_object.rds"))

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

UMAP_coord$cell_type = UMAP_coord$clusters
UMAP_coord[UMAP_coord$clusters == '1', 'cell_type'] = 'Unknown 1'
UMAP_coord[UMAP_coord$clusters == '3', 'cell_type'] = 'Unknown 2'
UMAP_coord[UMAP_coord$clusters == '5', 'cell_type'] = 'Early Germ Cells'
UMAP_coord[UMAP_coord$clusters == '6', 'cell_type'] = 'Middle Germ Cells 1'
UMAP_coord[UMAP_coord$clusters == '2', 'cell_type'] = 'Middle Germ Cells 2'
UMAP_coord[UMAP_coord$clusters == '4', 'cell_type'] = 'Late Germ Cells 1'


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
  scale_color_brewer(palette = 'Set2')
ggsave(filename = file.path(TARGET_dir, "celltypes.png"), plot = p, width = 8, height = 6)

norm_exp = monocle3::normalized_counts(cds)



UMAP_coord$bam = norm_exp['bam', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -bam), y=bam, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("bam normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_bam.png"), plot = p, width = 8, height = 6)

UMAP_coord$osk = norm_exp['osk', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -osk), y=osk, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("osk normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_osk.png"), plot = p, width = 8, height = 6)

UMAP_coord$osk = norm_exp['osk', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -osk), y=osk, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("osk normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_osk.png"), plot = p, width = 8, height = 6)

UMAP_coord$pgc = norm_exp['pgc', ] # quietscent 
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -pgc), y=pgc, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("pgc normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_pgc.png"), plot = p, width = 8, height = 6)

UMAP_coord$gcl = norm_exp['gcl', ] # quietscent 
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -gcl), y=gcl, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("gcl normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_gcl.png"), plot = p, width = 8, height = 6)

UMAP_coord$nos = norm_exp['nos', ] # quietscent 
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -nos), y=nos, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("nos normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_nos.png"), plot = p, width = 8, height = 6)

UMAP_coord$dpp = norm_exp['dpp', ]
p = ggplot(UMAP_coord, aes(x=reorder(cell_type, -dpp), y=dpp, fill = cell_type)) + 
  geom_violin() +
  guides(fill=guide_legend(title="")) +
  geom_boxplot(width=0.1) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set2') + 
  ylab("dpp normalized expression") + 
  xlab("cell type")
ggsave(filename = file.path(TARGET_dir, "violin_dpp.png"), plot = p, width = 8, height = 6)

###############################
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "4_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_late_GSEA_results.csv"))

focus_gsea = c("proteasomal protein catabolic process (GO:0010498)", 
               "proteasomal protein catabolic process (GO:0010498)", 
               "inorganic cation transmembrane transport (GO:0098662)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 6, 'Set2')[2]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw()
ggsave(filename = file.path(TARGET_dir, "late_Cells_GSEA_results.png"), plot = p, width = 10, height = 4)

# middle stage 2
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "2_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_middle_2_GSEA_results.csv"))

focus_gsea = c("proteasomal protein catabolic process (GO:0010498)", 
               "negative regulation of cellular macromolecule biosynthetic process (GO:2000113)", 
               "female gamete generation (GO:0007292)", 
               "negative regulation of translation (GO:0017148)", 
               "negative regulation of gene expression (GO:0010629)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 6, 'Set2')[4]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw()
ggsave(filename = file.path(TARGET_dir, "middle_2_Cells_GSEA_results.png"), plot = p, width = 10, height = 4)

# middle stage 1
GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "6_gsea_results.csv"), row.names = 1)
GSEA_results = GSEA_results[GSEA_results$padj < 0.05, ]
GSEA_results = GSEA_results[GSEA_results$NES > 0, ]
write.csv(GSEA_results, file = file.path(TARGET_dir, "sig_middle_1_GSEA_results.csv"))

focus_gsea = c("mRNA splicing, via spliceosome (GO:0000398)",
               "cytoskeleton-dependent cytokinesis (GO:0061640)", 
               "mitotic cytokinesis (GO:0000281)", 
               "cytoplasmic translation (GO:0002181)", 
               "DNA packaging (GO:0006323)")

sub_GSEA_results = GSEA_results[GSEA_results$pathway %in% focus_gsea, ]
sub_GSEA_results$log_pval = -log10(sub_GSEA_results$padj)

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 6, 'Set2')[3]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw()
ggsave(filename = file.path(TARGET_dir, "middle_1_Cells_GSEA_results.png"), plot = p, width = 10, height = 4)

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

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 6, 'Set2')[1]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw()
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

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 8, 'Set2')[5]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw()
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

p = ggplot(data=sub_GSEA_results, aes(x=reorder(pathway, log_pval), y=log_pval)) +
  geom_bar(stat="identity", fill = RColorBrewer::brewer.pal(n = 8, 'Set2')[6]) + coord_flip() + 
  xlab("GO Biological Processes") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("") +
  theme_bw()
ggsave(filename = file.path(TARGET_dir, "unknown_2_Cells_GSEA_results.png"), plot = p, width = 10, height = 4)

