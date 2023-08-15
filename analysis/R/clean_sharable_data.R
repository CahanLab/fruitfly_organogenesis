TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'clean_sharable_data')
dir.create(TARGET_dir, recursive = TRUE)

##### clean the stage 10-12 embryos #####
dir.create(file.path(TARGET_dir, 'stage_10_12'), recursive = TRUE)
object = readRDS(file.path('results', ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))
object@meta.data$Integrated_tentativeCellType = NULL
saveRDS(object, file = file.path(TARGET_dir, 'stage_10_12', 'stage10_12_seurat_object.rds'))

withr::with_dir(
  file.path(TARGET_dir, 'stage_10_12'), 
  {
    raw_matrix = object@assays$RNA@counts
    meta_data = object@meta.data
    write(colnames(raw_matrix), file = file.path("stage_10_12_cell_id_colnames.txt"))
    write(rownames(raw_matrix), file = file.path("stage_10_12_genes_rownames.txt"))
    Matrix::writeMM(raw_matrix, file.path("stage_10_12_raw_counts_sparse_matrix.txt"))
    write.csv(meta_data, file = 'stage_10_12_meta_data.csv')
  }
)

##### clean the stage 13-16 embryos #####
dir.create(file.path(TARGET_dir, 'stage_13_16'), recursive = TRUE)
object = readRDS(file.path('results', ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object4.rds"))
object@meta.data = object@meta.data[, -20:-48]
object@meta.data$RNA_snn_res.1.5 = NULL
saveRDS(object, file = file.path(TARGET_dir, 'stage_13_16', 'stage13_16_seurat_object.rds'))

withr::with_dir(
  file.path(TARGET_dir, 'stage_13_16'), 
  {
    raw_matrix = object@assays$RNA@counts
    meta_data = object@meta.data
    write(colnames(raw_matrix), file = file.path("stage_13_16_cell_id_colnames.txt"))
    write(rownames(raw_matrix), file = file.path("stage_13_16_genes_rownames.txt"))
    Matrix::writeMM(raw_matrix, file.path("stage_13_16_raw_counts_sparse_matrix.txt"))
    write.csv(meta_data, file = 'stage_13_16_meta_data.csv')
  }
)

##### clean data for SG monocle3 #####
cds = readRDS(file.path('results', ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland/monocle3_no_batch_correct_object.rds"))
cds@colData$subtypes = NA
cds@colData[clusters(cds) == 2, 'subtypes'] = 'Early Salivary Gland Cells'
cds@colData[clusters(cds) == 1, 'subtypes'] = 'Late Salivary Gland Cells'
dir.create(file.path(TARGET_dir, 'SG_specific'), recursive = TRUE)
saveRDS(object = cds, file = file.path(TARGET_dir, 'SG_specific', 'SG_monocle3_object.rds'))
##### clean data for trachea monocle3 #####
cds = readRDS(file.path('results', ANALYSIS_VERSION, "refined_wt_late_early_trachea/monocle3_no_batch_correct_object.rds"))
cds@colData$subtypes = NA
cds@colData[cds@colData$cell_type == 'Branching Trachea Cells', 'subtypes'] = 'Tracheal Tip Cells'
cds@colData[cds@colData$cell_type == 'Late Trachea Cells', 'subtypes'] = 'Late Tracheal Cells'
cds@colData[cds@colData$cell_type == 'Middle Trachea Cells', 'subtypes'] = 'Interm. Tracheal Cells'
cds@colData[cds@colData$cell_type == 'Early Trachea Cells', 'subtypes'] = 'Early Tracheal Cells'
cds@colData$cell_type = NULL  
dir.create(file.path(TARGET_dir, 'trachea_specific'), recursive = TRUE)
saveRDS(object = cds, file = file.path(TARGET_dir, 'trachea_specific', 'trachea_monocle3_object.rds'))
##### clean data for germ cells monocle3 #####
cds = readRDS(file.path('results', ANALYSIS_VERSION, "refined_wt_late_early_germ/monocle3_no_batch_correct_object.rds"))
cds@colData$subtypes = NA
cds@colData[clusters(cds) == '1', 'subtypes'] = 'Unknown 1'
cds@colData[clusters(cds) == '3', 'subtypes'] = 'Unknown 2'
cds@colData[clusters(cds) == '5', 'subtypes'] = 'Early Germ Cells'
cds@colData[clusters(cds) == '6', 'subtypes'] = 'Interm. Germ Cells 1'
cds@colData[clusters(cds) == '2', 'subtypes'] = 'Interm. Germ Cells 2'
cds@colData[clusters(cds) == '4', 'subtypes'] = 'Late Germ Cells'
dir.create(file.path(TARGET_dir, 'germ_specific'), recursive = TRUE)
saveRDS(object = cds, file = file.path(TARGET_dir, 'germ_specific', 'germ_monocle3_object.rds'))
