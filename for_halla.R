args = commandArgs(trailingOnly=TRUE)

# args[1]<-"/Users/feiwu/Downloads/ssgsea-results-scores.gct"
# args[2]<-"/Users/feiwu/Box/RNA-Seq/Dual-seq/samplesheet.csv"

setwd(dirname(args[1]))
setwd("../halla")
# read the original samplesheet
coldata0 <- read.csv(args[2], header = T, na.strings=c("","NA")) 
coldata<-coldata0[,1:2]
# read the ssGSEA score .gct
score<-read.table(args[1], row.names=1, sep="\t",header=T, skip = 2)
score<-score[coldata$sample_name]
# adjust covariance effect
score_adj <- limma::removeBatchEffect(score, coldata$group)
write.table(score_adj,"Host_score.txt",sep="\t",quote=F,col.names=NA)
