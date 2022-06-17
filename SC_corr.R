args = commandArgs(trailingOnly=TRUE) # Passing arguments to an R script from bash/shell command lines

filename<-basename(args[1]) # input; path to samplesheet_SC.csv
setwd(args[2]) # output directory
system(paste0("mkdir -p ",args[2],"/Correlations_MicrobiomeXhost"))

### Setup the Seurat objects ###
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
setwd(paste0(args[2],"/Correlations_MicrobiomeXhost"))

# Read samplesheet
samplesheet <- read.csv(args[1], header = T)
# Read location of host gene count matrix
host_matrix <- unique(samplesheet[2])

## function for uppercase the gene names
library(tibble)
library(dplyr)
Uppercase <- function(df){
  df<-add_column(df, NAME=toupper(row.names(df)), .before = 1)
  df<-df %>% group_by(NAME) %>% summarise_all(list(mean)) #average the duplicated rows
  df <- as.data.frame(df) # format tibble to dataframe
  row.names(df) <- df[,1] # rename the row name with uppercased gene names
  df <- df[,-1] # remove the "NAME" column
}

# function for QC
if (length(args)<7){
  SCT_m<-function(l){
    l <- PercentageFeatureSet(l,pattern = "^MT-", col.name = "percent.mt")
    l <- subset(l, subset= nFeature_RNA>200 & nFeature_RNA < 2*median(l@meta.data[,3]) & percent.mt < 10)
    l <- SCTransform(l, vars.to.regress = "percent.mt")
  }
} else {
  SCT_m<-function(l){
    l <- PercentageFeatureSet(l,pattern = "^MT-", col.name = "percent.mt")
    l <- subset(l, subset= nFeature_RNA > args[5] & nFeature_RNA < args[6] & percent.mt < args[7])
    l <- SCTransform(l, vars.to.regress = "percent.mt")
  }
}

if (nrow(host_matrix) > 1){
  SC.list <- list()
  for (h in host_matrix[,1]){
    if (args[3] == 1 || 3){ # 10x platform
      if (file.exists(paste0(h,"/matrix.mtx"))){
        print("Load the 10X datasets (matrix; folder)")
        host.df <- Read10X(h)
      } else if (length(list.files(h, pattern="\\.h5$"))>0){
        print("Load the 10X datasets (h5 format)")
        host.df <- Read10X_h5(paste0(h,"/",list.files(h, pattern="\\.h5$")))
      } else {print(paste0(h," : format of count matrix is not supported"))}
      # transfer to dataframe for reprocessing
      host.df <- as.data.frame(host.df)
      # clean the names of cell barcodes
      colnames(host.df)<-gsub("-1|.1","",colnames(host.df))
    } else if (args[3] == 2) { # Dropseq platform
      print ("Load the Dropseq datasets (.dge.txt)")
      host.df<-read.table(paste0(h,"/",list.files(h, pattern="\\.dge.txt$")),header = T, row.name=1)
      }
    # to uppercase the gene names
    host.df<-Uppercase(host.df)
    # Read microbiome SC data
    SC.micro<-read.delim(paste0(args[2],"/",basename(h),"_count_matrix.txt"), row.names = 1)
    # combine host and microbiome SC data
    metatranscripts<-rbind(host.df,SC.micro)
    # make dataframe to Seurat objects
    metatranscripts <- CreateSeuratObject(metatranscripts, project=basename(h),min.cells = 1, min.feature = 200)
    # add Seurat objects into a list
    SC.list <- append(SC.list, list(metatranscripts))
  } 
  # SCTransform for each dataset independently
  SC.list <- lapply(SC.list, SCT_m)
  features <- SelectIntegrationFeatures(object.list = SC.list, nfeatures = 3000)
  SC.list <- PrepSCTIntegration(object.list = SC.list, anchor.features = features)
  ### Perform integration ###
  SC.anchors <- FindIntegrationAnchors(object.list = SC.list, normalization.method = "SCT", anchor.features = features)
  SC.combined.sct <- IntegrateData(anchorset = SC.anchors, normalization.method = "SCT")
  
} else if (nrow(host_matrix) == 1){
  h <- host_matrix[,1]
  if (args[3] == 1 || 3){ # 10x platform
    if (file.exists(paste0(h,"/matrix.mtx"))){
      print("Load the 10X datasets (matrix; folder)")
      host.df <- Read10X(h)
    } else if (length(list.files(h, pattern="\\.h5$"))>0){
      print("Load the 10X datasets (h5 format)")
      host.df <- Read10X_h5(paste0(h,"/",list.files(h, pattern="\\.h5$")))
    } else {print(paste0(h," : format of count matrix is not supported"))}
    host.df <- as.data.frame(host.df)
    colnames(host.df)<-gsub("-1|.1","",colnames(host.df))
  } else if (args[3] == 2) { # Dropseq platform
    print ("Load the Dropseq datasets (.dge.txt)")
    host.df<-read.table(paste0(h,"/",list.files(h, pattern="\\.dge.txt$")),header = T, row.name=1)
  }
  host.df<-Uppercase(host.df)
  SC.micro<-read.delim(paste0(args[2],"/",basename(h),"_count_matrix.txt"), row.names = 1)
  metatranscripts<-rbind(host.df,SC.micro)
  metatranscripts <- CreateSeuratObject(metatranscripts, project=basename(h),min.cells = 1, min.feature = 200)
  metatranscripts <- SCT_m(metatranscripts)
  SC.combined.sct <- metatranscripts # unify the name to facilitate the downstream process
  features <- SC.combined.sct@assays[["SCT"]]@var.features
}

# dimensional reduction
SC.combined.sct <- RunPCA(SC.combined.sct, verbose = FALSE)
SC.combined.sct <- RunUMAP(SC.combined.sct, reduction = "pca", dims = 1:30)

# Find & label the clusters
SC.combined.sct <- FindNeighbors(SC.combined.sct, dims = 1:30)
SC.combined.sct <- FindClusters(SC.combined.sct, resolution = 0.5)

### Correlation test of microbiome and host gene among top 3,000 features ###
featured.microbiome <- subset(as.data.frame(SC.combined.sct@assays[["SCT"]]@data),row.names(SC.combined.sct@assays[["SCT"]]@data) %in% grep("taxid", features,value = T))
featured.host <- subset(as.data.frame(SC.combined.sct@assays[["SCT"]]@data),row.names(SC.combined.sct@assays[["SCT"]]@data) %in% grep("taxid", features,value = T, invert = T))
featured.microbiomeXhost <- cor(t(featured.microbiome),t(featured.host), method = "spearman")

# Parallel correlation test
library(parallel)
library(doParallel)
cores <- as.numeric(args[4])
options('mc.cores' = cores)
registerDoParallel(cores)
p.corr<-function(testmatrix1,testmatrix2) {
  x<-foreach (j = 1:ncol(testmatrix1),
              .combine = rbind,
              .multicombine = TRUE,
              .inorder = FALSE,
              .packages = c('data.table', 'doParallel')) %:%
    foreach (i = 1:ncol(testmatrix2),
             .combine = cbind,
             .multicombine = TRUE,
             .inorder = FALSE,
             .packages = c('data.table', 'doParallel')) %dopar% {
               a <- cor.test(testmatrix1[,j], testmatrix2[,i], method = "spearman") 
               a$p.value
             }
  colnames(x)<-colnames(testmatrix2)
  row.names(x)<-colnames(testmatrix1)
  return(x)
}

# Acquire the p-value of correlation test
p.mat<-p.corr(t(featured.host),t(featured.microbiome))

# Acquire the significant correlated microbiome and host genes
sig.corr<-data.frame()
for (i in 1:nrow(p.mat)){
  for (j in 1:ncol(p.mat)){
    if (p.mat[i,j]<0.05 && t(featured.microbiomeXhost)[i,j]>0.2){
      sig.corr <- rbind(sig.corr,data.frame(Microbiome=colnames(t(featured.microbiomeXhost))[j],
                                            Host_gene=rownames(t(featured.microbiomeXhost))[i],
                                            r=t(featured.microbiomeXhost)[i,j],
                                            p=p.mat[i,j]))
    }
  }
}
write.table(sig.corr, "MicrobiomeXhost_sigCorrelation.tsv", row.names = F, sep="\t",quote = F)

# Sort descending by correlations coefficient
sig.corr <- sig.corr[order(-sig.corr$r),]
# get the microbial names of top 10 highest correlations coefficient
sig.corr.top10<-sig.corr[1:10,]

# Show & Save the top 10 highest correlations
for (r in 1:nrow(sig.corr.top10)){
  FeaturePlot(object = SC.combined.sct, features = c(sig.corr.top10[r,2], sig.corr.top10[r,1]), blend = TRUE, order =T)
  ggsave(paste0(r,"_highest_cor.pdf"), width = 18, height = 6)
}

print("Significantly correlated microbiome and host genes are saved in /Correlations_MicrobiomeXhost/MicrobiomeXhost_sigCorrelation.tsv")
print("Top 10 highest correlations are saved in /Correlations_MicrobiomeXhost")
