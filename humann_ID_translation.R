args = commandArgs(trailingOnly=TRUE) # Passing arguments to an R script from bash/shell command lines
dir<-args[3]

translation<-function(tsv,terms){
  # import humann_kegg/go tables
  id<-read.table(tsv, sep="\t", comment.char = "", header = T, check.names = F)
  # adjust the column names
  names(id)[1] <- "Gene_Family"
  names(id)[-1] <- sub("_Abundance-RPKs","",names(id)[-1])
  # split KEGG/GO terms to a new column
  id$new <- sapply(regmatches(id$Gene_Family, gregexpr("^GO:[0-9]{7}|^K[0-9]{5}", id$Gene_Family)), paste, collapse = "")
  # pass ungrouped row name
  if ("" %in% id$new){id$new[match("",id$new)] <- id$Gene_Family[match("",id$new)]}
  # drop "ko:" row.names for KEGG terms table
  row.names(terms)<-gsub("^ko:","",row.names(terms))
  # left join by all.x=T, keep all data in "id" dataframe
  id1<-merge(id, terms, all.x=T, by.x="new", by.y="row.names")
  # pass unknown and ungrouped row name
  while (NA %in% id1$terms){id1$terms[match(NA,id1$terms)] <- id1$new[match(NA,id1$terms)]}
  row.names(id1)<-id1$terms
  id1<-subset(id1,select=-c(new,Gene_Family,terms))
  # output keeps first col rowname "Gene_Family"
  write.table(data.frame("Gene_Family"=rownames(id1),id1), sub(".tsv$","_translated.tsv",tsv), sep="\t",row.names = F)
}

## for KEGG
terms<-read.table(paste0(dir,"/koterms.txt"),sep="\t", quote="")
translation(args[1],terms)

## for GO
terms<-read.table(paste0(dir,"/goterms.txt"),sep="\t", quote="")
translation(args[2],terms)

