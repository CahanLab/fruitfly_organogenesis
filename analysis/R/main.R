# To reproduce our work, run this script. 
# We used Ubuntu 20.04 with an R package environment managed by renv. 
# We used R version 4.1.2.

# Make it work on each of our computers
try({setwd("~/Dropbox/fruitfly_organogenesis/analysis")}, silent = T)
try({setwd("~/Dropbox (CahanLab)/fruitfly_organogenesis/analysis")}, silent = T)
try({setwd("c:/Users/DAndrewLab/Desktop/fruitfly_organogenesis/analysis/")}, silent = T)
renv::activate()

ANALYSIS_VERSION = "v18" # Where to put output
metadata = read.table(header = T, text=
  "sample cellranger
   wt_rep1 2021-10-27_scRNA_10x_3prime 
   wt_rep3 2022-08-11_scRNA_10x_3prime_CellRanger6.1.2
   wt_early_rep1 2022-12-07_scRNA_10x_3prime_WT_early 
   wt_early_rep2 2022-12-13_scRNA_10x_3prime_WT_early2")

SAMPLE = "all"
source("R/set_up_environment.R")

for(i in 1:4){ 
  SAMPLE = metadata$sample[i]  
  CELLRANGER = metadata$cellranger[i] 
  source("R/set_up_environment.R") # Uses SAMPLE to set location of output
  source("R/filtering_and_qc.R") 
  source("R/refined_clustering.R")
  source("R/annotate_clusters_automated_BDGP.R") # perform a quick annotation using genes from BDGP 
}

# this is to integrate and cell type wild-type rep 1 and rep 3
source("R/integrate_wt13.R") 

# manually label the wt13 integrated data 
source("R/manual_label_wt13.R")

# plot out statistics for for wt13
source("R/plot_wt13_statistics.R")

# plot out UMAP and proportion for for wt13
source("R/plot_wt13_UMAP_proportion.R")

# cross study comparison for late wildtype
source("R/cross_study_comparison_wt13.R")

# plot out cross study comparison for late wildtype 
source("R/plot_cross_study_comparison_wt13.R")

# perform GSEA on the data 
source("R/wt13_GSEA.R")

# plot out the results for late wildtype 3 GSEA 
source("R/plot_wt13_GSEA.R")

# plot out the marker genes for rare cell types 
source("R/plot_selected_marker_genes_wt13.R")

# plot out all the marker genes for the cell types 
source("R/plot_all_marker_genes_wt13.R")

# this is to integrate and cell type wild-type early rep 1 and rep 2
source("R/integrate_early_wt12.R") 

# this is to manually label wildtype 1 and wildtype 2 early integrated data 
source("R/manual_label_early_wt12.R") 

# perform GSEA on the early wildtype data 
source("R/early_wt12_GSEA.R")

# plot out the results for late wildtype 3 GSEA 
source("R/plot_early_wt12_GSEA.R")

# plot out statistics for for early wt12
source("R/plot_early_wt12_statistics.R")

# plot out UMAP and proportion for for early wt12
source("R/plot_early_wt12_UMAP_proportion.R")

# cross study comparison for early wildtype
# will need to run this manually because you would have to run some python scripts 
source("R/cross_study_comparison_early_wt12.R")

#plot the cross study comparisons for early wild type 2
source("R/plot_cross_study_comparison_wt12.R")

# refined analysis of early + late salivary gland 
source("R/refined_wt13_early12_salivary_gland.R")

# make figures for the in depth analysis of salivary glands
source("R/plot_refined_wt13_early12_salivary_gland.R")

# refined analysis of early + late trachea cells 
source("R/refined_wt13_early12_trachea.R")

# make figures for the in depth analysis of trachea 
source("R/plot_refined_wt13_early12_trachea.R")

# refined analysis of early + late germ cells 
source("R/refined_wt13_early12_germ.R")

# make figures for the in depth analysis of germ cells 
source("R/plot_refined_wt13_early12_germ.R")

# match similarity between wild-type early samples and late samples 
source("R/match_early_clusters_to_late_clusters.R")

# compare cell cycle cells 
source("R/plot_compare_cell_cycle_wt13_early_12.R")

# plot out the quality metrics for all 4 batches 
source("R/plot_quality_metrics_4_batches.R")

# perform analysis of matrisome in late wildtype 
source("R/wt13_matrisome.R")

# perform analysis of matrisome in early wildtype 
source("R/early_wt12_matrisome.R")


######################################################
# down below are experimental scripts that I am trying out for testings. The results may or may not end up 
# going into the final manuscript 

# experimental script to compute the DE genes between early cluster and all late samples 
# this is to help verify some of the ambiguous early clusters 
source("R/early_cluster_vs_all_late.R")

# experimental script to compute the DE genes between all the clusters and identify gene enrichment 
# this is to help verify some of the ambiguous early clusters 
source("R/early_cluster_vs_others.R")

# test the germ cells 
source("R/germ_cells_exploration.R")

# match CNS subtypes with Seroka 
source("R/wt13_match_CNS_seroka.R")

# early wt12 celltyping compare with the celltyping in Calderon et al, 2022 DOI: 10.1126/science.abn5800
source("R/early_wt12_compare_Calderon.R")

# early wt12 celltyping compare with the celltyping in Seroka et al, 2022 DOI: https://doi.org/10.1016/j.ydbio.2022.05.018
source("R/early_wt12_compare_Seroka.R")

# late wt13 celltyping compare with the celltyping in Calderon et al, 2022 DOI: 10.1126/science.abn5800
source("R/wt13_compare_Calderon.R")

# late wt13 celltyping compare with the celltyping in Seroka et al, 2022 DOI: https://doi.org/10.1016/j.ydbio.2022.05.018
source("R/wt13_compare_Seroka.R")

