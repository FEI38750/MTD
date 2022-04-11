args = commandArgs(trailingOnly=TRUE)
library(tidyverse)

setwd(dirname(args[1]))
dir.create("../ssGSEA",recursive = T)
setwd("../ssGSEA")
# read the original samplesheet
coldata0 <- read.csv(args[2], header = T, na.strings=c("","NA")) 
coldata<-coldata0[,1:2]
# # extract the contrast information for reference
# coldata_vs<- coldata0[c("group1","group2")]
# coldata_vs<-coldata_vs[rowSums(is.na(coldata_vs)) == 0,] #remove the NA rows
# # read samplesheet as factors (as.is = F) for Deseq2 statistical analysis
# coldata_factor <- read.csv(args[2], header = T, as.is = F)
# coldata<-coldata_factor[,1:2]

gct<-read.csv(args[1],header = T)
gct1<-select(gct,"gene_name","hybrid_name",coldata[,1]) # select and reoder the columes
names(gct1)[1:2]<-c("NAME","Description")
gct1<-gct1[gct1["NAME"]!="-",]
# upcase the gene names
gct1$NAME<-toupper(gct1$NAME)

row1<-"#1.2"
row2<- c(nrow(gct1),length(coldata[,1]))
maxlength <- max(length(row1), length(row2), length(gct1))
length(row1) <- maxlength                   
length(row2) <- maxlength
gct2<-rbind(row1,row2,names(gct1),gct1)
gct2[1:2,][is.na(gct2[1:2,])]<-""

write.table(gct2,"host.gct",
            col.names = F,row.names = F,
            sep="\t",quote = F)
