library(singleCellNet)

TARGET_dir = file.path("results", ANALYSIS_VERSION, 'seroka_germ_classification')
dir.create(TARGET_dir, recursive = TRUE)

Seroka_object = readRDS("accessory_data/Doe_Drosophila_Embryo_Atlas/script/curated_embryo_Doe.rds")
our_monocle = readRDS("results/v18/refined_wt_late_early_germ/monocle3_no_batch_correct_object.rds")

UMAP_coord = our_monocle@int_colData$reducedDims$UMAP

colnames(UMAP_coord) = c("UMAP_1", "UMAP_2")
UMAP_coord = as.data.frame(UMAP_coord)
UMAP_coord$clusters = as.vector(monocle3::clusters(our_monocle))
UMAP_coord$pseudotime = as.vector(monocle3::pseudotime(our_monocle))
UMAP_coord$batch = as.vector(our_monocle@colData$batch)

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

train_exp = our_monocle@assays@data$counts
train_st = UMAP_coord

germ_Seroka_object = Seroka_object[, Seroka_object$cell_type == 'germline cells']
query_exp = germ_Seroka_object@assays$RNA@counts
query_st = germ_Seroka_object@meta.data

i_genes = intersect(rownames(train_exp), rownames(query_exp))

train_exp = train_exp[i_genes, ]
query_exp = query_exp[i_genes, ]

withr::with_dir(
  file.path(TARGET_dir), 
  { 
    write(colnames(train_exp), file = "raw_train_colnames.txt")
    write(rownames(train_exp), file = "raw_train_rownames.txt")
    Matrix::writeMM(train_exp, "raw_train_exp.txt")
    write.csv(train_st, file = 'raw_meta_tab.csv')
  }
)
withr::with_dir(
  file.path(TARGET_dir), 
  {
    write(colnames(query_exp), file = file.path("raw_query_colnames.txt"))
    write(rownames(query_exp), file = file.path("raw_query_rownames.txt"))
    Matrix::writeMM(query_exp, file.path("raw_query_exp.txt"))
  }
)

class_profile = read.csv(file.path(TARGET_dir, "SCN_classification.csv"), row.names = 1)
class_profile = class_profile[colnames(germ_Seroka_object), ]

germ_Seroka_object@meta.data = cbind(germ_Seroka_object@meta.data, class_profile)

DimPlot(germ_Seroka_object, group.by = 'SCN_class')
ggsave(filename = file.path(TARGET_dir, "SCN_class.png"))
FeaturePlot(germ_Seroka_object, features = 'Unknown.1')
ggsave(filename = file.path(TARGET_dir, "unknown_1.png"))

FeaturePlot(germ_Seroka_object, features = 'Unknown.2')
ggsave(filename = file.path(TARGET_dir, "unknown_2.png"))

FeaturePlot(germ_Seroka_object, features = 'Late.Germ.Cells')
ggsave(filename = file.path(TARGET_dir, "late_germ.png"))

FeaturePlot(germ_Seroka_object, features = 'Interm..Germ.Cells.2')
ggsave(filename = file.path(TARGET_dir, "interm2_germ.png"))

FeaturePlot(germ_Seroka_object, features = 'Interm..Germ.Cells.1')
ggsave(filename = file.path(TARGET_dir, "interm1_germ.png"))

FeaturePlot(germ_Seroka_object, features = 'Early.Germ.Cells')
ggsave(filename = file.path(TARGET_dir, "early_germ.png"))

germ_Seroka_object <- NormalizeData(germ_Seroka_object)
germ_Seroka_object <- FindVariableFeatures(germ_Seroka_object, selection.method = "vst", nfeatures = 2000)
germ_Seroka_object <- ScaleData(germ_Seroka_object)
germ_Seroka_object <- RunPCA(germ_Seroka_object, features = VariableFeatures(object = germ_Seroka_object))
germ_Seroka_object <- RunUMAP(germ_Seroka_object, dims = 1:10)

