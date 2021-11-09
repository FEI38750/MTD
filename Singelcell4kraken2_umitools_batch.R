#for run on server
args = commandArgs(trailingOnly=TRUE) # Passing arguments to an R script from command lines
library(reshape2)
files <- list.files(path=args[1], pattern="*.c.tsv$", full.names=TRUE, recursive=FALSE)
lapply(files, function(x) {
  t1<-read.table(x, sep="\t", quote="", head=T,comment.char = "")
  t1$cell<-gsub("b|'","",t1$cell)
  t2<-dcast(t1, gene ~ cell, value.var = "count" , fill=0)
  t3<-t2[,-1]
  rownames(t3)<-t2[,1] 
  write.table(data.frame("Name"=rownames(t3),t3),paste0(x,"_Count.txt"),sep="\t",row.names=FALSE)
})