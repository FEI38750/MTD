if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos = "http://cran.us.r-project.org")
BiocManager::install(version = "3.14")

BiocManager::install(c("biomaRt","DESeq2","tximeta","limma","phyloseq","glmGamPoi","cmapR","MAST",
                        "microbiome","ANCOMBC","Maaslin2","DO.db","clusterProfiler","enrichplot","pathview"))

if (!require("nloptr")) install.packages("nloptr",repos = "http://cran.us.r-project.org")

if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(tidyverse,ggrepel,colorspace,RColorBrewer,
                pheatmap,VennDiagram,doParallel,foreach,
                stringi,vegan,ggpubr,reshape2,sctransform,hdf5r,
                ggridges,ggnewscale,ggupset,Seurat)

if (!require("Seurat")) install.packages("Seurat",repos = "http://cran.us.r-project.org")

#included in tidyverse: ggplot2,tidyr,dplyr,stringr,tibble,