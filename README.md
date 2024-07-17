# Fruitfly Organogenesis Paper
This contains all the scripts and codes for the the paper [Organogenetic transcriptomes of the Drosophila embryo at single cell resolution](https://journals.biologists.com/dev/article/doi/10.1242/dev.202097/340234/Organogenetic-transcriptomes-of-the-Drosophila). 

## Processed Data 
The raw count matrices and processed data (Seurat object in .rds format) for stage 10-12 and stage 13-16 can be found on [GEO (GSE234602)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234602). The Seurat version used is 4.3.0. 
|Stage| Seurat object |
| --------- | --------------- |
| stage 10-12 | [download from GEO](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234602&format=file&file=GSE234602%5Fstage10%5F12%5Fseurat%5Fobject%2Erds%2Egz) |
| stage 13-16 | [download from GEO](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234602&format=file&file=GSE234602%5Fstage13%5F16%5Fseurat%5Fobject%2Erds%2Egz) |

The Monocle3 object (in .rds format) for the in-depth analysis of salivary gland, trachea and germ cells can be found below. The Monocle3 version is 1.3.1. 

|Cell Types | monocle3 object |
| --------- | --------------- |
| salivary gland | [download](https://cnobjects.s3.amazonaws.com/drosophila_2023/monocle3_objects/SG_monocle3_object.rds) |
| trachea | [download](https://cnobjects.s3.amazonaws.com/drosophila_2023/monocle3_objects/trachea_monocle3_object.rds) |
| germ cells | [download](https://cnobjects.s3.amazonaws.com/drosophila_2023/monocle3_objects/germ_monocle3_object.rds) |


## Navigating the analysis code 
All the analysis and figures presented in the manuscript are done on a 2020 Macbook Pro, M1 Chip with macOS Ventura 13.3.1 (a). To run the analysis code, you would need to download several raw data and accessory data files. 

1. Download and extract all the content from [here](https://cnobjects.s3.amazonaws.com/drosophila_2023/quantification.tar.gz). Place the folder `quantification` in the this directory. In short, these files are outputted by CellRanger. 

2. Download and extract all content from [here](https://cnobjects.s3.amazonaws.com/drosophila_2023/accessory_data.tar.gz) and place the folder called `accessory_data` in `analysis` directory creating a `/analysis/accessory_data` directory. In short, these files are accessory data (e.g matrisome genes, expression profiles from Seroka et al and Calderon et al) that help us perform the downstream analysis. 

3. After all the raw data and accessory data are downloaded and placed in the correct directory, please see [the README in analysis](/analysis/README.md) for more instructions on how to run the code. 
