library(reticulate)
library(sceasy)

use_condaenv('OneCC_dev')
loompy <- reticulate::import('loompy')

object = readRDS("results/v18/manual_annotation_wt13/manual_celltype_object4.rds")
sceasy::convertFormat(object, from="seurat", to="anndata",
                      outFile=file.path('results/v18/manual_annotation_wt13/', 'for_cellxgenes_object.h5ad'))

object = readRDS("results/v18/manual_annotation_early_wt12/manual_celltype_object1.rds")
sceasy::convertFormat(object, from="seurat", to="anndata",
                      outFile=file.path("results/v18/manual_annotation_early_wt12/", 'for_cellxgenes_object.h5ad'))
