# Fruitfly Organogenesis Paper
This contains all the scripts and codes for the single-cell fruit fly organogenesis paper (ref when published). 

## Processed Data 
The raw count matrices and processed data (Seurat object in .rds format) for stage 10-12 and stage 13-16 can be found on [GEO database]() (will release when the manuscript is accepted). 

The Monocle3 object (in .rds format) for the in-depth analysis of salivary gland, trachea and germ cells can be found below.


## Navigating the analysis code 
All the analysis and figures presented in the manuscript are done on a 2020 Macbook Pro, M1 Chip with macOS Ventura 13.3.1 (a). To run the analysis code, you would need to download several raw data and accessory data files. 

1. Create a folder named `quantification` in this directory and extract all the content from [here]() (will share when manuscript is accepted) into *quantification* directory. In short, these files contain the raw count matrices outputted by CellRanger. 

2. Create a folder named `accessory_data` in `analysis` directory and extract all content from [here]() (will share when manuscript is accepted) into `/analysis/accessory_data` directory. In short, these files are accessory data (e.g matrisome genes, expression profiles from Seroka et al and Calderon et al) that help us perform the downstream analysis. 

3. After all the raw data and accessory data are downloaded and placed in the correct directory, please see [the README in analysis](/analysis/README.md) for more instructions on how to run the code. 
