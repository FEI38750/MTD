if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos = "http://cran.us.r-project.org")
BiocManager::install(version = "3.12")

BiocManager::install(c("biomaRt","DESeq2","tximeta","limma","phyloseq","glmGamPoi","cmapR","MAST"))

if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(tidyverse,ggrepel,colorspace,RColorBrewer,
                pheatmap,VennDiagram,doParallel,foreach,
                stringi,vegan,ggpubr,reshape2,sctransform,hdf5r)

install.packages('Seurat',repos = "http://cran.us.r-project.org")

#included in tidyverse: ggplot2,tidyr,dplyr,stringr,tibble,