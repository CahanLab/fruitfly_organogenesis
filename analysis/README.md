# Where to find information
All the analysis and figures presented in the manuscript are done on a 2020 Macbook Pro, M1 Chip with macOS Ventura 13.3.1 (a). Down below are the softwares used for the main analysis

| Package     | Version     |
| ----------- | ----------- |
| R           | 4.1.2       |
| Seurat      | 4.3.0       |
| harmony     | 0.1.1       |
| monocle3    | 1.3.1       |
| presto      | 1.0.0       |
| fgsea       | 1.20.0      |

| Package     | Version     |
| ----------- | ----------- |
| Python      | 3.9.12      |
| pySCN       | commit [fb44024](https://github.com/pcahan1/PySingleCellNet/commit/fb440247f6f372367b8cca015c40b82f9f388573)         |

- Current code: `main.R` calls our current analysis scripts in the correct order. You can find all the R scripts in `R/` folder. I would recommend running the scripts sourced in `main.R` individually because `cross_study_comparison_early_wt12.R` and `cross_study_comparison_wt13.R` require you to run Python scripts somewhere in the middle. For the Python scripts that train on Seroka or Calderon data and apply to this study's you can find them [here](/python_scripts/pySCN_other/). You can find the Python scripts that train on our data and annotation and apply to Seroka or Calderon data [here](/python_scripts/pySCN_our/). Just place the scripts in the appropriate directory commented in `cross_study_comparison_early_wt12.R` and `cross_study_comparison_wt13.R`, and run them as instructed in the R scripts. 
- Current analysis: look for the highest version numbers in the `results` folder. The analysis version for the latest manuscript is V18. 
- There are also many R scripts in the Github that are not used in the manuscript but rather for general explorations of the data. 

