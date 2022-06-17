args = commandArgs(trailingOnly=TRUE) # Passing arguments to an R script from bash/shell command lines

filename<-basename(args[1]) # input; path to samplesheet_SC.csv
setwd(args[2]) # output directory
system(paste0("mkdir -p ",args[2],"/fastq"))
system(paste0("mkdir -p ",args[2],"/Cell_Barcode"))

samplesheet <- read.csv(args[1], header = T)
host_matrix <- unique(samplesheet[2])

# Concatenate the data from the same sample into one
for (i in 1:nrow(host_matrix)){
  u <- host_matrix[i,]
  for_cat <- samplesheet[,1][samplesheet[,2]==u]
  for_cat1 <- grep("1.fastq$",for_cat, value = T)
  for_cat2 <- grep("2.fastq$",for_cat, value = T)
  system(paste0("cat ",paste0(for_cat1, collapse =" ")," > ",u,"_1.fastq"))
  system(paste0("cat ",paste0(for_cat2, collapse =" ")," > ",u,"_2.fastq"))
  system(paste0("mv ",u,"_1.fastq ",args[2],"/fastq"))
  system(paste0("mv ",u,"_2.fastq ",args[2],"/fastq"))
}

if (args[3] == 1 || args[3] == 3){ # 10x platform
  # Make whitelist from host matrix (for h5 format)
  for (h in host_matrix[,1]){
    if (file.exists(paste0(h,"/matrix.mtx"))){
      print(paste0("this is matrix.mtx, whiltelist - barcodes.tsv should already exist in the folder: ",h))
      # Transfer whilelist to output/Cell_Barcode folder with name sampleName_Barcodes.tsv
      whitelist<-read.delim(paste0(h,"/",list.files(h, pattern="barcodes.tsv$")),header = F)
      whitelist<-gsub("-1","",whitelist[,1])
      write.table(whitelist,paste0(args[2],"/Cell_Barcode/",basename(h),"_Barcodes.tsv"),sep="\t",quote=F,row.names = F,col.names = F)
    } else if (length(list.files(h, pattern="\\.h5$"))>0){
      print(paste0("this is .h5 format, to make whiltelist - *_Barcodes.tsv in the folder: ",h))
      library(Seurat)
      h5_df <- Read10X_h5(paste0(h,"/",list.files(h, pattern="\\.h5$")))
      write.table(h5_df,paste0(h,"/",sub("\\.h5$",".tsv",list.files(h, pattern="\\.h5$"))),sep="\t",quote=F)
      # to extract the barcodes
      whitelist<-gsub("-1","",h5_df@Dimnames[[2]])
      write.table(whitelist,paste0(h,"/",sub("\\.h5$","_Barcodes.tsv",list.files(h, pattern="\\.h5$"))),sep="\t",quote=F,row.names = F,col.names = F)
      # copy whilelist to output/Cell_Barcode folder
      system(paste0("mv ",paste0(h,"/",list.files(h, pattern="_Barcodes.tsv$"))," ",args[2],"/Cell_Barcode/",basename(h),"_Barcodes.tsv"))
    } else {print(paste0(h," : Not supported to make whitlist, please check it manually"))}
  }
} else if (args[3] == 2) { # Dropseq platform
  for (h in host_matrix[,1]){
    dge.df<-read.table(paste0(h,"/",list.files(h, pattern="\\.dge.txt$")),header = T, row.name=1)
    write.table(colnames(dge.df),paste0(h,"/",sub("\\.dge.txt$","_Barcodes.tsv",list.files(h, pattern="\\.dge.txt$"))),sep="\t",quote=F,row.names = F,col.names = F)
    # copy whilelist to output/Cell_Barcode folder
    system(paste0("mv ",paste0(h,"/",list.files(h, pattern="_Barcodes.tsv$"))," ",args[2],"/Cell_Barcode/",basename(h),"_Barcodes.tsv"))
  }
}


