library(ggplot2)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'early_wt12_matrisome')

object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))

matrisome_df = read.csv("accessory_data/matrisome_data/drosophila_matrisome.csv")

interesting_cat_list = c("Basement Membrane", "Basement Membrane; Laminin", 'Insulin Family', 'Glue', "Other", "Lysyl Hydroxylase", "Prolyl 4-Hydroxylase")

# the below function as a modification from Seurat Dotplot https://github.com/satijalab/seurat/blob/HEAD/R/visualization.R
# only return the plotting DF 
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
 
big_plot_df = data.frame()
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
  #temp_plot_df$avg.exp.scaled = NULL
  temp_plot_df$matrisome_type = temp_cat 
  
  good_genes = filter_genes(temp_plot_df)
  temp_plot_df = temp_plot_df[temp_plot_df$features.plot %in% good_genes, ]
  
  big_plot_df = rbind(big_plot_df, temp_plot_df)
}
big_plot_df$log_exp = log1p(big_plot_df$avg.exp)
p <- ggplot(data = big_plot_df, mapping = aes_string(x = 'id', y = 'features.plot')) +
  geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp.scaled')) +
  #scale.func(range = c(0, 100), limits = c(scale.min, scale.max)) +
  guides(size = guide_legend(title = 'Percent Expressed')) +
  guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
  scale_colour_viridis_c() + 
  labs(
    x = 'Cell Types',
    y = 'Genes'
  ) + 
  theme_classic()  + 
  facet_grid(
    rows = vars(matrisome_type),
    scales = "free_y",
    space = "free_y",
    switch = "x"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.y = element_text(angle = 0), axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) +
  ggtitle("Stage 10-12 Embryos")
ggsave(filename = file.path(TARGET_dir, "plasmatocytes_genes.png"), plot = p, width = 10, height = 9)

##### This is to get the genes for Chitin cuticle genes #####
interesting_cat_list = c("Cuticle; Tweedle", "Cuticle", "Chitin-binding-domain-containing Proteins", "Cuticle; R&R Chitin-binding-domain-containing Proteins")
big_plot_df = data.frame()
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
  #temp_plot_df$avg.exp.scaled = NULL
  temp_plot_df$matrisome_type = temp_cat 
  
  good_genes = filter_genes(temp_plot_df)
  temp_plot_df = temp_plot_df[temp_plot_df$features.plot %in% good_genes, ]
  
  big_plot_df = rbind(big_plot_df, temp_plot_df)
}
big_plot_df$log_exp = log1p(big_plot_df$avg.exp)
p <- ggplot(data = big_plot_df, mapping = aes_string(x = 'id', y = 'features.plot')) +
  geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp.scaled')) +
  #scale.func(range = c(0, 100), limits = c(scale.min, scale.max)) +
  guides(size = guide_legend(title = 'Percent Expressed')) +
  guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
  scale_colour_viridis_c() + 
  labs(
    x = 'Cell Types',
    y = 'Genes'
  ) + 
  theme_classic()  + 
  facet_grid(
    rows = vars(matrisome_type),
    scales = "free_y",
    space = "free_y",
    switch = "x"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.y = element_text(angle = 0), axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) +
  ggtitle("Stage 10-12 Embryos")
ggsave(filename = file.path(TARGET_dir, "chitin_genes.png"), plot = p, width = 12, height = 19)

