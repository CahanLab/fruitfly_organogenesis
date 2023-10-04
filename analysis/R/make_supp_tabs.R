library(openxlsx)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'supplementary_tabs')
dir.create(TARGET_dir, recursive = TRUE)

##### output all the differentially expressed genes for each corresponding cell type #####
marker_gene_list = list()
marker_gene_tab = read.csv("results/v18/wt13_integrated/marker_genes.csv")
marker_gene_tab$X = NULL
marker_gene_tab$cell_types = NA

seurat_object = readRDS("results/v18/manual_annotation_wt13/manual_celltype_object4.rds")
seurat_meta = seurat_object@meta.data
for(temp_cluster in unique(marker_gene_tab$cluster)) {
  marker_gene_tab[marker_gene_tab$cluster == temp_cluster, 'cell_types'] = unique(seurat_meta[seurat_meta$seurat_clusters == temp_cluster, "manual_celltypes"]) 
}

marker_gene_list[['stage 13-16']] = marker_gene_tab
# write in the stage 10-12 
marker_gene_tab = read.csv("results/v18/early_wt12_integrated/marker_genes.csv")
marker_gene_tab$X = NULL
marker_gene_tab$cell_types = NA

seurat_object = readRDS("results/v18/manual_annotation_early_wt12/manual_celltype_object1.rds")
seurat_meta = seurat_object@meta.data
for(temp_cluster in unique(marker_gene_tab$cluster)) {
  marker_gene_tab[marker_gene_tab$cluster == temp_cluster, 'cell_types'] = unique(seurat_meta[seurat_meta$seurat_clusters == temp_cluster, "manual_celltypes"]) 
}
marker_gene_list[['stage 10-12']] = marker_gene_tab

write.xlsx(marker_gene_list, file = file.path(TARGET_dir, "Table_2.xlsx"))

##### output the SG, trachea and GC GSEA results #####
run_file_list = list()
run_file_list[['Stage13-16_SG']] = file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Salivary Gland/gsea_results_wt.csv")
run_file_list[['Stage13-16_Tr']] = file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Trachea/gsea_results_wt.csv")
run_file_list[['Stage13-16_GC']] = file.path("results", ANALYSIS_VERSION, "wt13_enrichment/Germ Cells/gsea_results_wt.csv")
run_file_list[['Stage10-12_SG']] = file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Salivary Gland/gsea_results_wt.csv")
run_file_list[['Stage10-12_Tr']] = file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Trachea/gsea_results_wt.csv")
run_file_list[['Stage10-12_GC']] = file.path("results", ANALYSIS_VERSION, "early_wt12_enrichment/Germ Cells/gsea_results_wt.csv")
table_list = list()
for(temp_name in names(run_file_list)) {
  gsea_results = read.csv(run_file_list[[temp_name]], row.names = 1)
  gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
  table_list[[temp_name]] = gsea_results
}
write.xlsx(table_list, file = file.path(TARGET_dir, "Table_4.xlsx"))

##### output GSEA results for SG subtypes #####
run_file_list = list()
run_file_list[['Early_SG']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "early_gsea_results.csv")
run_file_list[['Late_SG']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", "late_gsea_results.csv")
table_list = list()
for(temp_name in names(run_file_list)) {
  gsea_results = read.csv(run_file_list[[temp_name]])
  gsea_results$X = NULL
  gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
  table_list[[temp_name]] = gsea_results
}
write.xlsx(table_list, file = file.path(TARGET_dir, "SG_GSEA_subtypes.xlsx"))

##### output GSEA results for tracheal subtypes #####
run_file_list = list()
run_file_list[['Tip_Tr']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Branching Trachea Cells_gsea_results.csv")
run_file_list[['Early_Tr']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Early Trachea Cells_gsea_results.csv")
run_file_list[['Interm_Tr']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Middle Trachea Cells_gsea_results.csv")
run_file_list[['Late_Tr']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea", "Late Trachea Cells_gsea_results.csv")

table_list = list()
for(temp_name in names(run_file_list)) {
  gsea_results = read.csv(run_file_list[[temp_name]])
  gsea_results$X = NULL
  gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
  table_list[[temp_name]] = gsea_results
}
write.xlsx(table_list, file = file.path(TARGET_dir, "Tr_GSEA_subtypes.xlsx"))

##### output GSEA results for GC subtypes #####
run_file_list = list()
run_file_list[['Late_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "4_gsea_results.csv")
run_file_list[['Interm2_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "2_gsea_results.csv")
run_file_list[['Interm1_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "6_gsea_results.csv")
run_file_list[['Early_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "5_gsea_results.csv")
run_file_list[['Unknown1_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "1_gsea_results.csv")
run_file_list[['Unknown2_GC']] = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ/cluster_process_GSEA", "3_gsea_results.csv")

table_list = list()
for(temp_name in names(run_file_list)) {
  gsea_results = read.csv(run_file_list[[temp_name]])
  gsea_results$X = NULL
  gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
  table_list[[temp_name]] = gsea_results
}
write.xlsx(table_list, file = file.path(TARGET_dir, "GC_GSEA_subtypes.xlsx"))

##### make quality comparisons #####
output_statistics <- function(seurat_object) {
  return(c(median(seurat_object$nCount_RNA),
           median(seurat_object$nFeature_RNA), 
           mean(seurat_object$nCount_RNA),
           mean(seurat_object$nFeature_RNA)))
}

Seroka_object = readRDS("accessory_data/Doe_Drosophila_Embryo_Atlas/script/curated_embryo_Doe.rds")
Seroka_object = Seroka_object[, Seroka_object$dataset == 'stg12']
Calderon_object= readRDS("accessory_data/continuum_drosophila_embryonic_development_RNA/processed_data/continuum_exploration/10-12_celltyped.rds")
our_object = readRDS("results/v18/manual_annotation_early_wt12/manual_celltype_object1.rds")

stats_df = data.frame(row.names = c('median_nCount', 'median_nFeature', 'mean_nCount', 'mean_nFeature'))
stats_df$Seroka_stg12 = output_statistics(Seroka_object)
stats_df$Calderon_stg10_12 = output_statistics(Calderon_object)
stats_df$This_data_stg10_12 = output_statistics(our_object)

Seroka_object = readRDS("accessory_data/Doe_Drosophila_Embryo_Atlas/script/curated_embryo_Doe.rds")
Seroka_object = Seroka_object[, Seroka_object$dataset != 'stg12']
Calderon_object= readRDS("accessory_data/continuum_drosophila_embryonic_development_RNA/processed_data/continuum_exploration/14-16_celltyped.rds")
our_object = readRDS("results/v18/manual_annotation_wt13/manual_celltype_object4.rds")

stats_df$Seroka_stg14_16 = output_statistics(Seroka_object)
stats_df$Calderon_stg14_16 = output_statistics(Calderon_object)
stats_df$This_data_stg13_16 = output_statistics(our_object)

write.csv(stats_df, file = file.path(TARGET_dir, "quality_comparisons_tab.csv"))

##### make tables for genes for the different subtypes #####
# get the trachea 
DE_genes = read.csv("results/v18/refined_wt_late_early_trachea/rank_sum_test.csv", row.names = 1)
DE_genes[DE_genes$group == 'Branching Trachea Cells', 'group'] = 'Tracheal Tip Cells'
DE_genes[DE_genes$group == 'Late Trachea Cells', 'group'] = 'Late Tracheal Cells'
DE_genes[DE_genes$group == 'Middle Trachea Cells', 'group'] = 'Interm. Tracheal Cells'
DE_genes[DE_genes$group == 'Early Trachea Cells', 'group'] = 'Early Tracheal Cells'
write.csv(DE_genes, file.path(TARGET_dir, "trachea_subtype_genes.csv"))

# get the salivary 
DE_genes = read.csv("results/v18/refined_wt_late_early_salivary_gland/rank_sum_test.csv", row.names = 1)
DE_genes[DE_genes$group == 2, 'group'] = "Earlier Salivary Gland Cells"
DE_genes[DE_genes$group == 1, 'group'] = "Later Salivary Gland Cells"
write.csv(DE_genes, file.path(TARGET_dir, "salivary_subtype_genes.csv"))

# get the salivary 
DE_genes = read.csv("results/v18/refined_wt_late_early_germ/rank_sum_test.csv", row.names = 1)
DE_genes[DE_genes$group == '1', 'group'] = 'Unknown 1'
DE_genes[DE_genes$group == '3', 'group'] = 'Unknown 2'
DE_genes[DE_genes$group == '5', 'group'] = 'Early Germ Cells'
DE_genes[DE_genes$group == '6', 'group'] = 'Interm. Germ Cells 1'
DE_genes[DE_genes$group == '2', 'group'] = 'Interm. Germ Cells 2'
DE_genes[DE_genes$group == '4', 'group'] = 'Late Germ Cells'
write.csv(DE_genes, file.path(TARGET_dir, "germ_subtype_genes.csv"))

##### make table for matrisome (preparation) #####
modified_dotPlot_df <- function(
    object,
    assay = NULL,
    features,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = 'radius',
    scale.min = NA,
    scale.max = NA
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  
  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  return(data.plot)
}

filter_genes <- function(plot_df) {
  good_genes = c()
  for(temp_gene in unique(plot_df$features.plot)) {
    subset_plot_df = plot_df[plot_df$features.plot == temp_gene, ]
    if(max(subset_plot_df$pct.exp) > 1) {
      good_genes = c(good_genes, temp_gene)
    }
  }
  return(good_genes)
}

matrisome_df = read.csv("accessory_data/matrisome_data/drosophila_matrisome.csv")

##### to get the matrisome genes for stage10-12 #####
object = readRDS("results/v18/manual_annotation_early_wt12/manual_celltype_object1.rds")
gene_list = list()
gene_list[['collagen']] = c("vkg", "Col4a1")
gene_list[['Hydroxylase']] = c("Plod", "PH4alphaEFB")
gene_list[['laminins']] = c("LanA", "LanB1", "LanB2")
gene_list[['basement']] = c("Ppn", "Pxn", "SPARC", 'Tig')
gene_list[['phagocytic_receptor']] = c("drpr", 'crq', 'NimC4', 'Sr-CI', 'PGRP-LC')
gene_list[['Thioester-containing']] = c("NimB1", 'NimB3', "NimB4", "NimB5")
gene_list[['engulfment']] = c("SCAR", "Rac1", "Rac2", 'chic')
gene_list[['rabs']] = c("Rab5", "Rab7", "Lamp1")
gene_list[['cathepsin']] = c("cathD", "CtsB1")
gene_list[['insulin']] = c("Ilp4", "Ilp6")
gene_list[['antimicrobial']] = c("Drs", "Dro", "DptA", "DptB",'CecA1', "CecA2", 
                                 "Def", 'Mtk', 'BomS4', 'BomBc2', 'BomT1', 'BomS3', 
                                 'BomBc1', 'BomS2', 'BomT3', 'BomBc3', 'BomS6', 
                                 'BomS1', 'BomT2', 'BomS5')
genes_df = data.frame()
for(temp_cat in names(gene_list)) {
  temp_plot_df = modified_dotPlot_df(object, features = gene_list[[temp_cat]], group.by = 'manual_celltypes')
  temp_plot_df$matrisome_type = temp_cat 
  genes_df = rbind(genes_df, temp_plot_df)
}

interesting_cat_list = c("Cuticle; Tweedle", "Cuticle", "Chitin-binding-domain-containing Proteins", "Cuticle; R&R Chitin-binding-domain-containing Proteins")
for(temp_cat in interesting_cat_list) {
  if(temp_cat == 'Other') {
    matrisome_genes = c('Ppn', 'Pxn', 'SPARC', 'Tig')
  }
  else {
    temp_matrisome_df = matrisome_df[matrisome_df$Matrisome.Class...Protein.Family == temp_cat, ]
    matrisome_genes = temp_matrisome_df$Gene.Name
    matrisome_genes = intersect(matrisome_genes, rownames(object))
  }
  
  temp_plot_df = modified_dotPlot_df(object, features = matrisome_genes, group.by = 'manual_celltypes')
  temp_plot_df$matrisome_type = temp_cat 
  
  good_genes = filter_genes(temp_plot_df)
  temp_plot_df = temp_plot_df[temp_plot_df$features.plot %in% good_genes, ]
  
  genes_df = rbind(genes_df, temp_plot_df)
}

write.csv(genes_df, file = file.path(TARGET_dir, "stage10-12_matrisome.csv"))

##### to get the matrisome genes for stage13-16 #####
object = readRDS("results/v18/manual_annotation_wt13/manual_celltype_object4.rds")
gene_list = list()
gene_list[['collagen']] = c("vkg", "Col4a1")
gene_list[['Hydroxylase']] = c("Plod", "PH4alphaEFB")
gene_list[['laminins']] = c("LanA", "LanB1", "LanB2")
gene_list[['basement']] = c("Ppn", "Pxn", "SPARC", 'Tig')
gene_list[['phagocytic_receptor']] = c("drpr", 'crq', 'NimC4', 'Sr-CI', 'PGRP-LC')
gene_list[['Thioester-containing']] = c("NimB1", 'NimB3', "NimB4", "NimB5")
gene_list[['engulfment']] = c("SCAR", "Rac1", "Rac2", 'chic')
gene_list[['rabs']] = c("Rab5", "Rab7", "Lamp1")
gene_list[['cathepsin']] = c("cathD", "CtsB1")
gene_list[['insulin']] = c("Ilp4", "Ilp6")
gene_list[['antimicrobial']] = c("Drs", "Dro", "DptA", "DptB",'CecA1', "CecA2", 
                                 "Def", 'Mtk', 'BomS4', 'BomBc2', 'BomT1', 'BomS3', 
                                 'BomBc1', 'BomS2', 'BomT3', 'BomBc3', 'BomS6', 
                                 'BomS1', 'BomT2', 'BomS5')
genes_df = data.frame()
for(temp_cat in names(gene_list)) {
  temp_plot_df = modified_dotPlot_df(object, features = gene_list[[temp_cat]], group.by = 'manual_celltypes')
  temp_plot_df$matrisome_type = temp_cat 
  genes_df = rbind(genes_df, temp_plot_df)
}

interesting_cat_list = c("Cuticle; Tweedle", "Cuticle", "Chitin-binding-domain-containing Proteins", "Cuticle; R&R Chitin-binding-domain-containing Proteins")
for(temp_cat in interesting_cat_list) {
  if(temp_cat == 'Other') {
    matrisome_genes = c('Ppn', 'Pxn', 'SPARC', 'Tig')
  }
  else {
    temp_matrisome_df = matrisome_df[matrisome_df$Matrisome.Class...Protein.Family == temp_cat, ]
    matrisome_genes = temp_matrisome_df$Gene.Name
    matrisome_genes = intersect(matrisome_genes, rownames(object))
  }
  
  temp_plot_df = modified_dotPlot_df(object, features = matrisome_genes, group.by = 'manual_celltypes')
  temp_plot_df$matrisome_type = temp_cat 
  
  good_genes = filter_genes(temp_plot_df)
  temp_plot_df = temp_plot_df[temp_plot_df$features.plot %in% good_genes, ]
  
  genes_df = rbind(genes_df, temp_plot_df)
}

write.csv(genes_df, file = file.path(TARGET_dir, "stage13-16_matrisome.csv"))

##### identify BDGP cell types for selected marker genes  #####

# load in the dictionary to convert the gene names to flybase name 

stage10_12_markers = read.csv("accessory_data/manual_marker_genes_v18/manual_marker_stage_10-12.csv")
stage10_12_markers$flybase_id = NULL
stage10_12_markers$BDGP_cell_type = NULL

#' @title Auto assign celltype from BDGP reference by accessing the website 
#' @description
#' Automatically assigns a cell type for a cluster based on DE genes and BDGP in-situs reference by accessing the 
#' BDGP website directly. This is just a programmatic way to copy all the cell types annotation from BDGP systematically without 
#' me going BDGP and manually copy the cell types.  
#' @param stage_name the stage that we are interested in 
#' @param gene_name flybase ID 
#' @return a vector of cell types that have expression in gene of interest 
BDGP_robot <- function(stage_name = 'stage13-16', gene_name = 'FBgn0003254') { 
  # set up a sleep function before we access the website to avoid overloading the server 
  Sys.sleep(5)
  url = paste0('https://insitu.fruitfly.org/cgi-bin/ex/report.pl?ftype=1&ftext=', gene_name)
  fly_data <- url %>%
    rvest::read_html() %>%
    rvest::html_nodes(xpath='/html/body/div[1]/div[2]/div[1]/table[2]') %>%
    rvest::html_nodes("tbody") %>% rvest::html_children()
  
  if(length(fly_data) == 0) { 
    return(vector())
  }
  
  for(temp_index in seq(1, length(fly_data))) { 
    cur_stage = fly_data[[temp_index]] %>% 
      rvest::html_nodes("th") %>% 
      rvest::html_text()
    if(cur_stage == stage_name) { 
      item = fly_data[[temp_index]] %>% 
        rvest::html_nodes("td") %>% 
        rvest::html_text()
      item_string = paste(item, collapse = '')
      item_string = stringr::str_remove_all(item_string, "\t")
      cell_types = stringr::str_split(item_string, "\n")[[1]]
      cell_types = cell_types[cell_types != ""]
      cell_types = trimws(cell_types)
    }
  }
  return(unique(cell_types))
}

for(gene_interest in stage10_12_markers$gene.abbr.) { 
  flybase_interest = as.character(gene_converter[gene_interest])
  stage10_12_markers[stage10_12_markers$gene.abbr. == gene_interest, 'flybase_id'] = flybase_interest
  BDGP_cellTypes = BDGP_robot(stage_name = 'stage11-12', gene_name = gene_interest)
  
  print(gene_interest)

  stage10_12_markers[stage10_12_markers$gene.abbr. == gene_interest, 'BDGP_cell_type'] = paste(BDGP_cellTypes, collapse = ";")
}

write.csv(stage10_12_markers, file = file.path(TARGET_dir, "stage10_12_markers_bdgp.csv"))

##### make the tables for stage 13-16 markers 
stage13_16_markers = read.csv("accessory_data/manual_marker_genes_v18/manual_marker_stage_13-16.csv")
stage13_16_markers$flybase_id = NULL
stage13_16_markers$BDGP_cell_type = NULL

for(gene_interest in stage13_16_markers$gene.abbr.) { 
  flybase_interest = as.character(gene_converter[gene_interest])
  stage13_16_markers[stage13_16_markers$gene.abbr. == gene_interest, 'flybase_id'] = flybase_interest
  BDGP_cellTypes = BDGP_robot(stage_name = 'stage13-16', gene_name = gene_interest)
  
  print(gene_interest)
  
  stage13_16_markers[stage13_16_markers$gene.abbr. == gene_interest, 'BDGP_cell_type'] = paste(BDGP_cellTypes, collapse = ";")
}

write.csv(stage13_16_markers, file = file.path(TARGET_dir, "stage13_16_markers_bdgp.csv"))
