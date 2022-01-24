args = commandArgs(trailingOnly=TRUE)

setwd(dirname(args[1]))
setwd("../halla")
# read the original samplesheet
coldata0 <- read.csv(args[2], header = T, na.strings=c("","NA")) 
coldata<-coldata0[,1:2]
# update the coldata if metadata is provided
if (length(args) == 3){
  coldata <- read.csv(args[3], header = T, as.is = F)
  coldata[]<-lapply(coldata, factor)
}
# read the ssGSEA score .gct
score<-read.table(args[1], row.names=1, sep="\t",header=T, skip = 2)
score<-score[coldata$sample_name]
# adjust covariance effect
score_adj <- limma::removeBatchEffect(score, coldata$group)
if (length(args) == 3){
  coldata.n<-coldata
  coldata.n[]<-lapply(coldata.n, as.numeric)
  normtrans_adj <- limma::removeBatchEffect(normtrans, covariates=coldata.n[,2:ncol(coldata.n)])
}
write.table(score_adj,"Host_score.txt",sep="\t",quote=F,col.names=NA)
