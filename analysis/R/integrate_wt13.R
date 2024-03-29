# integrate stage 13-16 wild type samples 
# we are integrating replicate 1 and replicate 3 because replicate 2 has quality issues 
# in the manuscript, replicate 3 is renamed as replicate 2 to avoid confusions 
SAMPLETYPE = "wt"
TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt13_integrated")
dir.create(TARGET_dir)

rep1 = readRDS(file.path("results", ANALYSIS_VERSION, paste0(SAMPLETYPE, "_rep1"), "BDGP_automated_annotation_object.rds"))
rep1@meta.data$batch = 'rep_1'

rep3 = readRDS(file.path("results", ANALYSIS_VERSION, paste0(SAMPLETYPE, "_rep3"), "BDGP_automated_annotation_object.rds"))
rep3@meta.data$batch = 'rep_3'

merged_seurat = merge(rep1, rep3)

raw_object = Seurat::CreateSeuratObject(merged_seurat@assays$RNA@counts, project = 'merged_samples')
raw_object@meta.data = merged_seurat@meta.data

saveRDS(raw_object, file = file.path(TARGET_dir, paste0(SAMPLETYPE, "_raw_merged.rds")))

object = Seurat::NormalizeData(raw_object)
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

genes_in_object = rownames(object@assays$RNA@meta.features)

mitochondrially_encoded_genes = read.table("accessory_data/mitochondrially_encoded_genes.tsv", header = F)[[1]]
ribosomal_rna_genes           = read.table("accessory_data/ribosomal_rna_genes.tsv", header = F)[[1]]
ribosomal_protein_genes       = grep("^rpl|^rps", ignore.case = T, value = T, genes_in_object)
object[["mitochondrial_transcript_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( convert_fbgn_to_symbol( mitochondrially_encoded_genes ), ) %>%
  colSums
object[["ribosomal_rna_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( convert_fbgn_to_symbol( ribosomal_rna_genes ), ) %>%
  colSums
object[["ribosomal_protein_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( ribosomal_protein_genes, ) %>%
  colSums
object[["mitochondrial_transcript_norm_expression"]] = object[["mitochondrial_transcript_total_expression"]] / object$nCount_RNA
object[["ribosomal_rna_norm_expression"           ]] = object[["ribosomal_rna_total_expression"           ]] / object$nCount_RNA
object[["ribosomal_protein_norm_expression"       ]] = object[["ribosomal_protein_total_expression"       ]] / object$nCount_RNA
object[["log10_nCount_RNA"]] = object[["nCount_RNA"]] %>% log10

object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
object %<>% ScaleData(features = VariableFeatures(object = object))
object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)
pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(object)), pdim=length(VariableFeatures(object)), var=1)
singular_values = slot(Reductions(object, slot = "pca"), "stdev")
is_significant = singular_values^2 > pc_cutoff
num_pc = sum(is_significant)

# after getting the significant PCs, rerun PCA 
object %<>% RunPCA(features = VariableFeatures(object = object), npcs = num_pc)
object %<>% FindNeighbors(dims = 1:num_pc)
object %<>% RunUMAP(dim = 1:num_pc)
object %<>% FindClusters()
pdf(file = file.path(TARGET_dir, "naive_integration.pdf"))
DimPlot(object, group.by = "batch")
dev.off()

# harmony integration 
object %<>% harmony::RunHarmony("batch")
object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc)
object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
object %<>% FindClusters()
pdf(file = file.path(TARGET_dir, "harmony_integration.pdf"))
DimPlot(object, group.by = "batch")
dev.off()

saveRDS(object, file = file.path(TARGET_dir, 'object.rds'))

# cell type the clusters using BDGP 
object = readRDS(file.path(TARGET_dir, 'object.rds'))

# extract top 20 cluster markers
withr::with_dir(
  file.path(TARGET_dir), 
  {
    if('top_marker_genes_1-20.csv' %in% list.files() == FALSE) { 
      differentially_expressed = FindAllMarkers(object, test.use = "bimod")
      write.csv(differentially_expressed, "marker_genes.csv")
      differentially_expressed %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(n=20, wt = avg_log2FC) %>%
        dplyr::arrange(cluster) %>%
        write.csv("top_marker_genes_1-20.csv")      
    } else { 
      differentially_expressed = read.csv("top_marker_genes_1-20.csv", row.names = 1)
    }
  }
)

# extract top 50 cluster markers
withr::with_dir(
  file.path(TARGET_dir), 
  {
    if('top_marker_genes_1-50.csv' %in% list.files() == FALSE) { 
      differentially_expressed = read.csv("marker_genes.csv", row.names = 1)
      differentially_expressed %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(n=50, wt = avg_log2FC) %>%
        dplyr::arrange(cluster) %>%
        write.csv("top_marker_genes_1-50.csv")      
    } 
  }
)

# load in the necessary code to automatically cell type the clusters via BDGP database 
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
gene_converter = as.vector(gene_metadata$flybase_id)
names(gene_converter) = as.vector(gene_metadata$symbol)

BDGP_database_ct = read.csv("accessory_data/BDGP_marker_genes/insitu_annot.csv", header = FALSE)
colnames(BDGP_database_ct) = c("gene_name1", "gene_name2", 'fly_base_id', 'stage', 'cell_type')

BDGP_database_image = read.csv("accessory_data/BDGP_marker_genes/insitu_images.csv", header = FALSE)
BDGP_database_image = BDGP_database_image[, seq(1, 7)]

withr::with_dir(
  file.path(TARGET_dir), 
  {
    top_marker_genes = read.csv("top_marker_genes_1-20.csv")
    top_marker_genes$cell_type = NULL
    top_marker_genes$flybase_id = NULL
    for(temp_index in rownames(top_marker_genes)) { 
      gene_interest = top_marker_genes[temp_index, 'gene']
      flybase_interest = as.character(gene_converter[gene_interest])
      top_marker_genes[temp_index, 'flybase_id'] = flybase_interest
      BDGP_cellTypes = BDGP_ct_assign(BDGP_database_ct, gene_name = flybase_interest)
      
      print(gene_interest)
      print(paste(BDGP_cellTypes, collapse = ";"))
      
      top_marker_genes[temp_index, 'cell_type'] = paste(BDGP_cellTypes, collapse = ";")
    }
    
    write.csv(top_marker_genes, "top_marker_genes_1-20_cellTypes.csv")
    
    cluster_vector = vector()
    cell_type_vector = vector()
    num_genes_vector = vector()
    gene_scores_vector = vector()
    
    for(temp_cluster in unique(top_marker_genes$cluster)) {
      temp_top_marker = top_marker_genes[top_marker_genes$cluster == temp_cluster, ]
      num_genes_list = list()
      score_list = list()
      
      for(temp_index in rownames(temp_top_marker)) { 
        log2Fold = temp_top_marker[temp_index, "avg_log2FC"]
        possible_cellTypes = stringr::str_split(temp_top_marker[temp_index, "cell_type"], ";")[[1]]
        if(possible_cellTypes[1] != "") { 
          for(possible_cellType in possible_cellTypes) { 
            if(possible_cellType %in% names(num_genes_list)) {
              num_genes_list[[possible_cellType]] = num_genes_list[[possible_cellType]] + 1
              score_list[[possible_cellType]] = score_list[[possible_cellType]] + log2Fold
            }
            else { 
              num_genes_list[[possible_cellType]] = 1
              score_list[[possible_cellType]] = log2Fold
            }
          }
        }
      }
      
      num_genes_list = unlist(num_genes_list)
      score_list = unlist(score_list)
      
      # remove the no staining 
      num_genes_list = num_genes_list[names(num_genes_list) != 'no staining']
      score_list = score_list[names(score_list) != 'no staining']
      
      # remove the no ubiquitous 
      num_genes_list = num_genes_list[names(num_genes_list) != 'ubiquitous']
      score_list = score_list[names(score_list) != 'ubiquitous']
      
      # remove the genes that are None 
      num_genes_list = num_genes_list[names(num_genes_list) != 'None']
      score_list = score_list[names(score_list) != 'None']
      
      num_genes_list = sort(num_genes_list, decreasing = TRUE)
      score_list = sort(score_list, decreasing = TRUE)
      
      cluster_vector = c(cluster_vector, temp_cluster)
      if(num_genes_list[names(score_list)[1]] < 2) { 
        cell_type_vector = c(cell_type_vector, 'Unknown')
        num_genes_vector = c(num_genes_vector, 0)
        gene_scores_vector = c(gene_scores_vector, 0)
      } else{ 
        if(score_list[names(score_list)[1]] < 5) {
          cell_type_vector = c(cell_type_vector, 'Unknown')
          num_genes_vector = c(num_genes_vector, num_genes_list[names(score_list)[1]])
          gene_scores_vector = c(gene_scores_vector, score_list[names(score_list)[1]])
        }
        else{
          cell_type_vector = c(cell_type_vector, names(score_list)[1])
          num_genes_vector = c(num_genes_vector, num_genes_list[names(score_list)[1]])
          gene_scores_vector = c(gene_scores_vector, score_list[names(score_list)[1]])
        }
      }
    }
    
    clusterAnnotation = data.frame(cluster = cluster_vector, 
                                   n_supporting_genes = num_genes_vector, 
                                   cellType_score = gene_scores_vector, 
                                   annotation = cell_type_vector)
    write.csv(clusterAnnotation, "tentativeCellTypes_BDGP.csv")
  }
)

withr::with_dir(
  file.path(TARGET_dir), 
  {
    top_marker_genes = read.csv('top_marker_genes_1-20_cellTypes.csv', row.names = 1)
    top_marker_genes$Images = NULL
    for(temp_index in rownames(top_marker_genes)) { 
      gene_interest = top_marker_genes[temp_index, 'gene']
      flybase_interest = as.character(gene_converter[gene_interest])
      top_marker_genes[temp_index, 'flybase_id'] = flybase_interest
      image_links = BDGP_image_assign(BDGP_database_image, gene_name = flybase_interest)
      
      print(gene_interest)
      
      top_marker_genes[temp_index, 'Images'] = image_links
    }
    top_marker_genes = subset(top_marker_genes, select = c('avg_log2FC', 'p_val_adj', 'cluster', 'gene', 'flybase_id', 'cell_type', 'Images'))
    top_marker_genes_df = DT::datatable(top_marker_genes, escape = FALSE, fillContainer = TRUE, filter = 'top', style = 'auto', width = '100px', height = '200px')
    doc <- htmltools::tagList(
      htmltools::div(top_marker_genes_df, style = "width=auto;height=auto")
    )
    doc[[1]]$children[[1]]$width = 'auto'
    doc[[1]]$children[[1]]$height = '800px'
    htmltools::save_html(html = doc, file = "widgets.html")
  }
)

# get tentative cell types 
object@meta.data$tentativeCellType = NULL
for(temp_cluster in unique(clusterAnnotation$cluster)) { 
  object@meta.data[object@meta.data$seurat_clusters == temp_cluster, 'Integrated_tentativeCellType'] = clusterAnnotation[clusterAnnotation$cluster == temp_cluster, 'annotation']
}

# remove the embryonic status 
object@meta.data$Integrated_tentativeCellType = stringr::str_remove(object@meta.data$Integrated_tentativeCellType, "embryonic ")
object@meta.data$Integrated_tentativeCellType = stringr::str_remove(object@meta.data$Integrated_tentativeCellType, "embryonic/larval ")

DimPlot(object, group.by = "Integrated_tentativeCellType", label = T, label.size = 5) +
  coord_fixed() + 
  ggtitle("Tentative cell type")
ggsave(file.path(TARGET_dir, "tentativeCellTypes_BDGP.png" ), width = 20, height = 7 )

DimPlot(object, group.by = "seurat_clusters", label = T, label.size = 5) +
  coord_fixed() + 
  ggtitle("seurat_clusters")
ggsave(file.path(TARGET_dir, "seurat_clusters.png" ), width = 20, height = 7 )

saveRDS(object, file = file.path(TARGET_dir, 'BDGP_automated_annotation_object.rds'))


