
TARGET_dir = file.path("results", ANALYSIS_VERSION, "x_chromosome_germ_cells")
dir.create(TARGET_dir)

##### load in all the data #####
monocle3_obj = readRDS("results/v18/refined_wt_late_early_germ/monocle3_no_batch_correct_object.rds")

monocle3_obj@colData[clusters(monocle3_obj) == 1, 'subtypes'] = 'Unknown 1'
monocle3_obj@colData[clusters(monocle3_obj) == 3, 'subtypes'] = 'Unknown 2'
monocle3_obj@colData[clusters(monocle3_obj) == 5, 'subtypes'] = 'Early Germ Cells'
monocle3_obj@colData[clusters(monocle3_obj) == 6, 'subtypes'] = 'Interm. Germ Cells 1'
monocle3_obj@colData[clusters(monocle3_obj) == 2, 'subtypes'] = 'Interm. Germ Cells 2'
monocle3_obj@colData[clusters(monocle3_obj) == 4, 'subtypes'] = 'Late Germ Cells'

wt_late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object4.rds"))
sub_wt_late_object = subset(wt_late_object, subset = manual_celltypes == "Germ Cells")
sub_wt_late_object$experimental_condition = 'late'

sub_wt_late_object <- FindNeighbors(sub_wt_late_object, dims = 1:66)
sub_wt_late_object <- FindClusters(sub_wt_late_object, resolution = 0.15)
DimPlot(sub_wt_late_object)
ggsave(filename = file.path(TARGET_dir, "new_cluster_label_late.png"))

wt_early_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))
sub_wt_early_object = subset(wt_early_object, subset = manual_celltypes == 'Germ Cells')
sub_wt_early_object$experimental_condition = 'early'
DimPlot(sub_wt_early_object)
ggsave(filename = file.path(TARGET_dir, "new_cluster_label_early.png"))


##### perform comparison of marker genes #####

sheet_index_list = c(1, 2, 3, 4) 
for(sheet_index in sheet_index_list) {
  if(sheet_index == 1) {
    type = "sexed_female"
  } else if(sheet_index == 2) {
    type = "sexed_male"
  } else if(sheet_index == 3) {
    type = "unsexed_female"
  } else if(sheet_index == 4) {
    type = "unsexed_male"
  }
  sex_marker_genes = readxl::read_excel("accessory_data/sex_determining_genes/SupplementalTable_S4.xlsx", sheet = sheet_index)
  DotPlot(sub_wt_early_object, features = sex_marker_genes$gene_short_name, scale = FALSE) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(file.path(TARGET_dir, paste0(type, "_early.png")), width = 9)
  
  FeaturePlot(sub_wt_early_object, features = sex_marker_genes$gene_short_name)
  ggsave(file.path(TARGET_dir, paste0(type, "_UMAP_early.png")), width = 20, height = 20)
  
  DotPlot(sub_wt_late_object, features = sex_marker_genes$gene_short_name, scale = FALSE)+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(file.path(TARGET_dir, paste0(type, "_late.png")), width = 9)
  
  FeaturePlot(sub_wt_late_object, features = sex_marker_genes$gene_short_name)
  ggsave(file.path(TARGET_dir, paste0(type, "_UMAP_late.png")), width = 20, height = 20)
  
  monocle3::plot_genes_by_group(monocle3_obj, markers = sex_marker_genes$gene_short_name, group_cells_by = 'subtypes') + coord_flip()
  ggsave(file.path(TARGET_dir, paste0(type, "_GC.png")), width = 9)
  monocle3::plot_cells(monocle3_obj, genes = as.vector(sex_marker_genes$gene_short_name), cell_size = 1, show_trajectory_graph = FALSE)
  ggsave(file.path(TARGET_dir, paste0(type, "_UMAP_GC.png")), width = 20, height = 20)
}

##### get X chromosome genes #####
genes_chrom = read.csv("../quantification/reference_genome_info/dmel-all-r6.33.gtf", header = FALSE, sep = '\t')

get_x_ratio <- function(seurat_obj, genes_chrom) {
  chrom_labels = c("3L", "3R", '2L', "2R", '4', "X", 'Y')
  gene_chrom_exp = matrix(data = NA, nrow = length(unique(seurat_obj$seurat_clusters)), ncol = length(chrom_labels))
  rownames(gene_chrom_exp) = unique(seurat_obj$seurat_clusters)
  colnames(gene_chrom_exp) = chrom_labels
  for(chrom_label in chrom_labels) {
    sub_genes_chrom = genes_chrom[genes_chrom$V1 == chrom_label, ]
    genes_symbol = stringr::str_split(sub_genes_chrom$V9, ";", simplify = TRUE)
    genes_symbol = as.data.frame(genes_symbol)
    genes_symbol = stringr::str_remove_all(genes_symbol$V2, " gene_symbol ")
    interest_genes = unique(genes_symbol)
    
    for(cluster in unique(seurat_obj$seurat_clusters)) {
      temp_object = subset(seurat_obj, subset = seurat_clusters == cluster)
      temp_exp = temp_object@assays$RNA@data
      temp_exp = temp_exp[intersect(rownames(temp_exp), interest_genes), ]
      gene_chrom_exp[cluster, chrom_label] = mean(apply(temp_exp, FUN = mean, MARGIN = 1))
      
    }
  }
  return(gene_chrom_exp)
}

gene_ratio = get_x_ratio(sub_wt_late_object, genes_chrom)
gene_ratio['1', ] / gene_ratio['0', ]

gene_ratio = get_x_ratio(sub_wt_early_object, genes_chrom)
gene_ratio['32', ] / gene_ratio['24', ]
