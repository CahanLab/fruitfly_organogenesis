# To reproduce our work, run this script. 
# We used MacOS 13.3.1 (a) with an R package environment managed by renv. 
# We used R version 4.1.2.

# Make it work on each of our computers
try({setwd("~/Dropbox/fruitfly_organogenesis/analysis")}, silent = T)
try({setwd("~/Dropbox (CahanLab)/fruitfly_organogenesis/analysis")}, silent = T)
try({setwd("c:/Users/DAndrewLab/Desktop/fruitfly_organogenesis/analysis/")}, silent = T)
renv::activate()

ANALYSIS_VERSION = "v18" # Where to put output

# wildtype rep 3 is presented as wildtype stage 13-16 rep 2 in the manuscript 
metadata = read.table(header = T, text=
  "sample cellranger
   wt_rep1 2021-10-27_scRNA_10x_3prime 
   wt_rep3 2022-08-11_scRNA_10x_3prime_CellRanger6.1.2
   wt_early_rep1 2022-12-07_scRNA_10x_3prime_WT_early 
   wt_early_rep2 2022-12-13_scRNA_10x_3prime_WT_early2")

SAMPLE = "all"
source("R/set_up_environment.R")

# do a quick automated preprocessing of the samples individually 
for(i in 1:4){ 
  SAMPLE = metadata$sample[i]  
  CELLRANGER = metadata$cellranger[i] 
  source("R/set_up_environment.R") # Uses SAMPLE to set location of output
  source("R/filtering_and_qc.R") 
  source("R/refined_clustering.R")
  source("R/annotate_clusters_automated_BDGP.R") # perform a quick annotation using genes from BDGP 
}

##### the main ordering of scripts to reproduce the figures/analysis in the manuscript #####
# this is to integrate and cell type wild-type rep 1 and rep 3 (in the manuscript rep 3 is labelled as rep 2) 
source("R/integrate_wt13.R") 

# manually label the wt13 integrated data 
source("R/manual_label_wt13.R")

# plot out UMAP and proportion plots for for wt13
# Fig 1b, 1c
# Supp Fig 1a, 1b
source("R/plot_wt13_UMAP_proportion.R")

# perform GSEA on the data 
source("R/wt13_GSEA.R")

# plot out the results for late wildtype 13 GSEA 
# Fig. 2D 
source("R/plot_wt13_GSEA.R")

# plot out the marker genes for rare cell types 
# Supp Fig. 17 - 28 
source("R/plot_selected_marker_genes_wt13.R")

# this is to integrate and cell type wild-type early rep 1 and rep 2 
source("R/integrate_early_wt12.R") 

# this is to manually label wildtype 1 and wildtype 2 early integrated data 
source("R/manual_label_early_wt12.R") 

# perform GSEA on the early wildtype data 
source("R/early_wt12_GSEA.R")

# plot out the results for early wildtype GSEA 
# Fig 2D
source("R/plot_early_wt12_GSEA.R")

# plot out UMAP and proportion for for early wt12 
# Fig 1D, 1E
# Supp Fig 1C, 1D
source("R/plot_early_wt12_UMAP_proportion.R")

# plot out the marker genes for rare cell types 
# Supp Fig 13-16
source("R/plot_selected_marker_genes_early_wt12.R")

# refined analysis of early + late salivary gland 
source("R/refined_wt13_early12_salivary_gland.R")

# make figures for the in depth analysis of salivary glands 
# Fig 3
# Supp Fig 2
source("R/plot_refined_wt13_early12_salivary_gland.R")

# refined analysis of early + late trachea cells 
source("R/refined_wt13_early12_trachea.R")

# make figures for the in depth analysis of trachea 
source("R/plot_refined_wt13_early12_trachea.R")

# refined analysis of early + late germ cells 
source("R/refined_wt13_early12_germ.R")

# make figures for the in depth analysis of germ cells 
# Fig 5 
# Supp Fig 5, 6
source("R/plot_refined_wt13_early12_germ.R")

# plot out the sequencing depth for all 4 batches 
source("R/plot_quality_metrics_4_batches.R")

# perform analysis of matrisome in late wildtype 
source("R/wt13_matrisome.R")

# perform analysis of matrisome in early wildtype 
source("R/early_wt12_matrisome.R")

# plot out matrisome genes 
# Fig 6A 
# Supp Fig 7A, 9
source("R/plot_wt12_matrisome.R")

# plot out matrisome genes 
# Fig 6B
# Supp Fig 7B, 10
source("R/plot_wt13_matrisome.R")

# correlation cell lines bulk expression with single cell clusters 
# Fig 7A
# Supp Fig 11 
source("R/assess_celllines.R")

# make the supplementary tables and tables shown in the manuscript
source("R/make_supp_tabs.R")

# identify the expressed genes in Seroka et al 
# Fig 7D
source("R/seroka_sg_genes.R")

# clean up the scRNA-seq data for sharing
source("R/clean_sharable_data.R")

##### scripts to run line by line #####
# these scripts require the user to run python scripts somewhere in the middle
# I could have streamline this better, but I am lazy. 
# cross study comparison for stage 13-16 wildtype 
# run this script individually 
source("R/cross_study_comparison_wt13.R")

# plot out cross study comparison for late wildtype 
# Fig. 7c
# Supp Fig 12b, 30, 31
source("R/plot_cross_study_comparison_wt13.R")

# cross study comparison for early wildtype 
# will need to run this manually because you would have to run some python scripts in the middle 
source("R/cross_study_comparison_early_wt12.R")

#plot the cross study comparisons for early wild type 2 
# Fig 7C
# Supp Fig 12A, 29
source("R/plot_cross_study_comparison_wt12.R")


##### unused scripts #####
# here are the scripts that were used for general exploration, 
# but did not make it in the main manuscript. 
# These scripts are also semi-interesting enough that I don't want to delete them. 

# verify whether or not the unknown clusters found in our germ cell 
# populations are also present in seroka et al data
source("R/check_seroka_germ.R")

# convert the Seurat objects (R) into anndata (python)
# we used Seurat to anndata mainly to help us explore the data interactively 
# on CellXGenes 
source("R/convert_to_h5ad.R")

# to explore the x chromosome gene expressions in our germ cell populations 
# and also explore the sex determination genes found in 
# https://genome.cshlp.org/content/31/6/1011.long
source("R/x_chromosome_exp.R")

# find the genes expressed in SG classified cells from Seroka et al data
source("R/seroka_sg_genes.R")

# make the plots for the quality of the 4 batches 
source("R/plot_quality_metrics_4_batches.R")

# match cell typing results of stage 10-12 with cell typing results of stage 13-16 
source("match_early_clusters_to_late_clusters.R")