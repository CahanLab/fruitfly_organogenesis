
TARGET_dir = file.path("results", ANALYSIS_VERSION, "x_chromosome_germ_cells")

monocle3_obj = readRDS("results/v18/refined_wt_late_early_germ/monocle3_no_batch_correct_object.rds")
female_markers = c("fs(1)h", 'Dsp1', 'Ubi-p5E', 'Pits', 'Smr', 'Imp', 'CG11699', 'N', 'CG32767', 'CG15891', 'mRpL14', 'rudhira', 'mRpL22', 'CG2865', 'Tis11', 'CG3224')
male_markers = c("vig2", 'Hsp27', 'Hsp23')
GC_markers = c("ovo", 'Trf2', 'bip2', 'Tif-IA', 'Cnot4', 'Not3', 'mamo')

monocle3::plot_genes_by_group(monocle3_obj, markers = female_markers)
monocle3::plot_genes_by_group(monocle3_obj, markers = male_markers)

wt_late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object4.rds"))
sub_wt_late_object = subset(wt_late_object, subset = manual_celltypes == "Germ Cells")
sub_wt_late_object$experimental_condition = 'late'

wt_early_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12/manual_celltype_object1.rds"))
sub_wt_early_object = subset(wt_early_object, subset = manual_celltypes == 'Germ Cells')
sub_wt_early_object$experimental_condition = 'early'
DotPlot(sub_wt_early_object, features = female_markers, scale = FALSE)
DotPlot(sub_wt_early_object, features = male_markers, scale = FALSE)

DotPlot(sub_wt_early_object, features = 'Ccdc114', scale = FALSE)

x_chrom_markers = c("be", 'CG1572', 'CG12773', 'CG16721', 'CG34332', 'Mnr', 'OtopLc', 'Rip11', 'Tob', 'zyd')
DotPlot(sub_wt_early_object, features = x_chrom_markers, scale = FALSE)

male_markers = c("dup", "side-V", 'tHMG1')
male_markers = c("Ccdc114", "CG3251", "CG12477", "CG13577", "CG18605", 'CG42857', 'Droj2', 'EIF4E6', 'Hsp83', "La", "mus301", 'NANS', "Rrp5", "Sf3b3", 'side-V', "srl")

DotPlot(sub_wt_early_object, features = male_markers, scale = FALSE)
