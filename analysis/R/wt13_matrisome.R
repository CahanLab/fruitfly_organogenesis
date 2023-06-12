# a systematic screening of matrisome genes 
# if you somehow arrived at this script, feel free to explore other genes in the matrisome 
matrisome_df = read.csv("accessory_data/matrisome_data/drosophila_matrisome.csv")

TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt13_matrisome")
dir.create(TARGET_dir)

object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13", "manual_celltype_object4.rds"))

# curate data 
make_output_df <- function(object, ECM_type, ct_col = 'manual_celltypes') {
  output_df = data.frame('cell_type' = unique(object@meta.data[, ct_col]), 
                         'ECM_type' = ECM_type)
  output_df$module_score = NA
  rownames(output_df) = output_df$cell_type
  
  for(celltype in unique(object@meta.data[, ct_col])) {
    module_score = mean(object@meta.data[object@meta.data[, ct_col] == celltype, 'expression_module_score1'])
    output_df[celltype, 'module_score'] = module_score
  }
  
  output_df$scaled_module_score = (output_df$module_score - mean(output_df$module_score)) / sd(output_df$module_score)
  return(output_df)
}


big_df = data.frame()
for(ECM_type in unique(matrisome_df$Matrisome.Category)) {
  sub_matrisome = matrisome_df[matrisome_df$Matrisome.Category == ECM_type, ]
  sub_genes = sub_matrisome$Gene.Name
  sub_genes = intersect(sub_genes, rownames(object))
  object = AddModuleScore(object, 
                          features = list(sub_genes),
                          name='expression_module_score')
  temp_df = make_output_df(object, ECM_type)
  big_df = rbind(big_df, temp_df)
}

p = ggplot(big_df, aes(cell_type, ECM_type, fill= scaled_module_score)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_viridis(discrete=FALSE) + 
  xlab("Cell Type (Stage 13-16)") + 
  ylab("Matrisome Type") + 
  labs(fill='Scaled Gene Module Score') 
ggsave(file.path(TARGET_dir, "Matrisome_category_score.png"), plot = p, height = 4, width = 10)

##### convert some secreted data as response to virus #####
matrisome_df[matrisome_df$Gene.Ontology._Cellular.component == 'defense response to virus [GO:0051607]', "Matrisome.Class...Protein.Family"] = 'Denfense Response to Virus'
matrisome_df = matrisome_df[matrisome_df$Matrisome.Class...Protein.Family != '', ]
for(temp_class in unique(matrisome_df$Matrisome.Category)) {
  sub_matrisome_df = matrisome_df[matrisome_df$Matrisome.Category == temp_class, ]
  big_df = data.frame()
  for(ECM_type in unique(sub_matrisome_df$Matrisome.Class...Protein.Family)) {
    sub_matrisome = matrisome_df[matrisome_df$Matrisome.Class...Protein.Family == ECM_type, ]
    sub_genes = sub_matrisome$Gene.Name
    sub_genes = intersect(sub_genes, rownames(object))
    if(length(sub_genes) == 0) {
      next
    }
    object = AddModuleScore(object, 
                            features = list(sub_genes),
                            name='expression_module_score')
    temp_df = make_output_df(object, ECM_type)
    big_df = rbind(big_df, temp_df)
  }
  
  p = ggplot(big_df, aes(cell_type, ECM_type, fill= scaled_module_score)) + 
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    scale_fill_viridis(discrete=FALSE) + 
    xlab("Cell Type (Stage 13-16)") + 
    ylab("Matrisome Type") + 
    labs(fill='Scaled Gene Module Score') + 
    ggtitle(temp_class)
  ggsave(file.path(TARGET_dir, paste0(temp_class, "_protein_family_score.png")), plot = p, height = 4, width = 10)
}

##### plot out the individual genes in the each category #####
matrisome_df = read.csv("accessory_data/matrisome_data/drosophila_matrisome.csv")
for(category in unique(matrisome_df$Matrisome.Category)) {
  dir.create(file.path(TARGET_dir, category))
  sub_matrisome_df = matrisome_df[matrisome_df$Matrisome.Category == category, ]
  for(temp_class in unique(sub_matrisome_df$Matrisome.Class...Protein.Family)) {
    sub_matrisome_df_2 = sub_matrisome_df[sub_matrisome_df$Matrisome.Class...Protein.Family == temp_class, ]
    if(temp_class == '') {
      title = 'unlabelled'
    } else {
      tite = temp_class
    }
    igenes = intersect(sub_matrisome_df_2$Gene.Name, rownames(object))
    if(length(igenes) == 0) {
      next
    }
    p = DotPlot(object, features = igenes, group.by = 'manual_celltypes') + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      ggtitle(tite) + coord_flip()
    ggsave(file.path(TARGET_dir, category,paste0(tite, "_dot_plot.png")), plot = p, height = 7, width = 10)
  }
}