object=readRDS(file.path(RESULTS, "object.Rdata"))
Idents(object)=object[["seurat_clusters"]][[1]]

# extract cluster markers
withr::with_dir(
  file.path(RESULTS), 
  {
    if('top_marker_genes_1-40.csv' %in% list.files() == FALSE) { 
      differentially_expressed = FindAllMarkers(object, test.use = "bimod")
      write.csv(differentially_expressed, "marker_genes.csv")
      differentially_expressed %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(n=20, wt = avg_log2FC) %>%
        dplyr::arrange(cluster) %>%
        write.csv("top_marker_genes_1-40.csv")      
    } else { 
      differentially_expressed = read.csv("top_marker_genes_1-40.csv", row.names = 1)
    }
  }
)

# load in the dictionary to convert the gene names to flybase name 
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

# based on the marker genes in BDGP and the expression of differentially expressed genes 
# assign a tentative cell type 
withr::with_dir(
  file.path(RESULTS), 
  {
    top_marker_genes = read.csv("top_marker_genes_1-40.csv")
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
      
      # remove the ubiquitous 
      num_genes_list = num_genes_list[names(num_genes_list) != 'ubiquitous']
      score_list = score_list[names(score_list) != 'ubiquitous']
      
      # remove the ubiquitous 
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

# to add the image links for each gene 
# construct a web-app of top DE genes for each cluster of cells. 
withr::with_dir(
  file.path(RESULTS), 
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

object@meta.data$tentativeCellType = NULL
for(temp_cluster in unique(clusterAnnotation$cluster)) { 
  object@meta.data[object@meta.data$seurat_clusters == temp_cluster, 'tentativeCellType'] = clusterAnnotation[clusterAnnotation$cluster == temp_cluster, 'annotation']
}

# remove 'embryonic'  
object@meta.data$tentativeCellType = stringr::str_remove(object@meta.data$tentativeCellType, "embryonic ")
object@meta.data$tentativeCellType = stringr::str_remove(object@meta.data$tentativeCellType, "embryonic/larval ")

DimPlot(object, group.by = "tentativeCellType", label = T, label.size = 5) +
  coord_fixed() + 
  ggtitle("Tentative cell type")
ggsave(file.path( RESULTS, "tentativeCellTypes_BDGP.png" ), width = 20, height = 7 )

saveRDS(object, file = file.path( RESULTS, 'BDGP_automated_annotation_object.rds'))
