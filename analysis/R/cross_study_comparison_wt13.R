library(Seurat)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "cross_study_comparison_wt13")
dir.create(TARGET_dir)

our_data = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_wt13/manual_celltype_object4.rds'))