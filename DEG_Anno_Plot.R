args = commandArgs(trailingOnly=TRUE) # Passing arguments to an R script from bash/shell command lines
# make folder DEG to store outputs
setwd(dirname(args[1]))
#system("mkdir -p DEG") # make new a folder DEG in current working directory

library(DESeq2)
library(tibble)
# Read & preprocess the input file cts before run deseq2 analysis
filename<-basename(args[1])
if (filename %in% c("humann_genefamilies_Abundance_go_translated.tsv","humann_genefamilies_Abundance_kegg_translated.tsv")){
  # read the data file
  #translated file has additional quote symbol "\""
  cts<-read.table(args[1],
                  row.names=1, sep="\t",header=T)
  #round the values before as.integer
  cts<-as.data.frame(lapply(lapply(cts,round),as.integer),row.names = row.names(cts))
  # make a folder for outputs
  if (filename == "humann_genefamilies_Abundance_go_translated.tsv"){
    dir.create("Nonhost_hmn_DEG/GO",recursive = T)
    setwd("Nonhost_hmn_DEG/GO")
  }
  if (filename == "humann_genefamilies_Abundance_kegg_translated.tsv"){
    dir.create("Nonhost_hmn_DEG/KEGG",recursive = T)
    setwd("Nonhost_hmn_DEG/KEGG")
  }
}
if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
  # read the data file
  #bracken file (eg. bracken_species_all) without quote symbol "\""; mark empty quote as quote=""
  cts<-read.table(args[1],
                  row.names=1, sep="\t",header=T, quote="")
  # Extract columns with count numbers
  cts<-cts[,grepl("*_num", names(cts))]
  # match column name to sample name
  colnames(cts)<-gsub("^Report_|\\.species.bracken_num$|\\.phylum.bracken_num$|.genus.bracken_num$","",colnames(cts))
  # decontamination step
  # read a list of contamination organisms
  # conta_ls<-read.table("~/Dual-seq/conta_ls.txt",sep="\t",colClasses = c("character","NULL"))
  # conta_ls<-conta_ls[conta_ls != "",]
  # remove the entries matching the contamination organisms
  # for (c in conta_ls){
  #   cts<-cts[!grepl(c,row.names(cts)),]
  # }
  # make a folder for outputs
  dir.create("Nonhost_DEG", recursive = T)
  setwd("Nonhost_DEG")
  # save the decontaminated result for later reference
  # write.table(cts,paste0(filename,"_decontaminated"),sep="\t", quote = F, col.names = NA)
}
if (filename == "host_counts.txt"){
  # read the data file
  #featureCounts file (eg. host_counts.txt) without quote symbol "\""; mark empty quote as quote=""
  cts<-read.table(args[1],
                  row.names=1, sep="\t",header=T, quote="")
  # drop first 5 columns with information other than counts
  cts<-cts[,-c(1:5)]
  # drop rows of zero count
  cts<-cts[rowSums(cts[-1])>0,]
  # make a folder for outputs
  dir.create("Host_DEG", recursive = T)
  setwd("Host_DEG")
}

# read the original samplesheet
coldata0 <- read.csv(args[2], header = T, na.strings=c("","NA"))
# extract the contrast information for reference
coldata_vs<- coldata0[c("group1","group2")]
coldata_vs<-coldata_vs[rowSums(is.na(coldata_vs)) == 0,] #remove the NA rows
# read samplesheet as factors (as.is = F) for Deseq2 statistical analysis
coldata_factor <- read.csv(args[2], header = T, as.is = F)
coldata<-coldata_factor[,1:2]

if (filename == "host_counts.txt"){
  # make cts(count matrix) has consistent order with samplesheet
  cts<-cts[coldata$sample_name]
  # load the datastructure to DESeq
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design= ~ group)
}

if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all",
                    "humann_genefamilies_Abundance_go_translated.tsv","humann_genefamilies_Abundance_kegg_translated.tsv")){
  if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
    files_h <- list.files(path=paste0(dirname(args[1]),"/temp"), pattern="^Report_host_.*\\.txt$", full.names=TRUE, recursive=FALSE)
  }
  if (filename %in% c("humann_genefamilies_Abundance_go_translated.tsv","humann_genefamilies_Abundance_kegg_translated.tsv")){
    humann_f <- gsub("/hmn_genefamily_abundance_files","",dirname(args[1]))
    files_h <- list.files(path=paste0(humann_f,"/temp"), pattern="^Report_host_.*\\.txt$", full.names=TRUE, recursive=FALSE)
  }
  # to get a list for host
  transcriptome_size<-c() #generate a empty list
  for (i in files_h){
    t<-read.table(i,sep="\t", quote= "")
    total_reads<-t[1,2] + t[2,2] #get total reads abundance of a sample
    fn<-gsub("^Report_host_|\\.txt$","",basename(i)) #grep the sample name
    transcriptome_size<-c(transcriptome_size,setNames(total_reads,fn)) #add sample name with value to the list lh
  }
  transcriptome_size <- log2(transcriptome_size)-mean(log2(transcriptome_size))
  coldata$order<-1:nrow(coldata)
  coldata<-merge(coldata,as.data.frame(transcriptome_size), by.x="sample_name",by.y="row.names")
  coldata<-coldata[order(coldata$order), ]
  coldata<-subset(coldata, select = -c(order))
  # make cts(count matrix) has consistent order with samplesheet
  cts<-cts[coldata$sample_name]
  # load the datastructure to DESeq
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design= ~ group + transcriptome_size)
}

# perform the DESeq analysis
dds <- DESeq(dds)  


# Data transformation for visualization (normalization included)
#rld<-rlog(dds,blind=F) # regularized log transformation (log2 based)
if (dim(results(dds))[1]  < 1000 || min(colSums(cts !=0)) < 1000){
  vsd<-varianceStabilizingTransformation(dds,blind=F) # vatiance stabilizing transformation
} else {
  vsd<-vst(dds,blind=F)
}
normtrans<-assay(vsd)

# save normalized & transformed data for visualization
write.csv(normtrans,file=paste0(sub(".tsv$|.txt$","",filename),"_normalized_transformed.csv"))

# save normalized (untransformed) data for reference
norm<-counts(dds,normalized=T)
write.csv(norm,file=paste0(sub(".tsv$|.txt$","",filename),"_normalized.csv"))

# merge and add suffixes; normalized and normalized&transformed
merge.nt<-merge(norm,normtrans,by="row.names", suffixes=c(".norm",".normtrans"))
if (filename == "host_counts.txt"){
  # add annotations; try up to 5 times if biomaRt server not response
  counter<-0
  while (exists("gene_ID")==F & counter<=5){
    library("biomaRt")
    if (args[3]=="9544"){
      ensembl=useMart("ensembl",dataset="mmulatta_gene_ensembl")
    }
    if (args[3]=="9606"){
      ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
    }
    if (args[3]=="10090"){
      ensembl=useMart("ensembl",dataset="mmusculus_gene_ensembl")
    }
    names(merge.nt)[1]<-"GeneID"
    genes <- merge.nt$GeneID
    gene_ID <- getBM(filters="ensembl_gene_id",
                     attributes=c("external_gene_name","ensembl_gene_id",
                                  "chromosome_name","start_position","end_position",
                                  "strand","gene_biotype","description"),
                     values=genes,mart=ensembl)
    names(gene_ID)[names(gene_ID)=="external_gene_name"]<-"gene_name" #rename the first column
    counter<-counter+1
  }
}
# #to add entrezid separately; caution: may bring a issue of "duplicate" ENTREZID!!
# library(clusterProfiler)
# library(org.Mmu.eg.db)
# # may bring a issue of "duplicate" ENTREZID!!
# ENTREZID = bitr(gene_ID$external_gene_name,
#                 fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mmu.eg.db") #get entrezid with acceptable id duplicate
# gene_ID <- merge(gene_ID,ENTREZID,by.x="gene_name",by.y="SYMBOL",all.x = T)

# function for host
# function of contrast between groups, merge to annotation and count tables, and save to .csv files
comparison<-function(dds,coldata_vs,filename,gene_ID,cts,merge.nt){
  for (i in 1:nrow(coldata_vs)){
  group1<-coldata_vs$group1[i]
  group2<-coldata_vs$group2[i]
  comparison <- results(dds, contrast=c("group", group1, group2))
  comparison_f<-as.data.frame(comparison)
  comparison_f<-comparison_f[order(comparison_f$pvalue),]
  system(paste0("mkdir -p"," ",group1,"_vs_",group2)) #for output file structure
  setwd(paste0(group1,"_vs_",group2))
  write.csv(comparison_f,file=paste0(sub(".tsv$|.txt$","",filename),"_",group1,"_vs_",group2,".csv"))
  setwd("../")
  colnames(comparison_f)<-paste0(colnames(comparison_f),".",group1,"_vs_",group2)
  gene_ID<-merge(gene_ID,comparison_f,by.x="ensembl_gene_id", by.y="row.names") #merge each comparison to annotation
  }
  Anno.merge.all<-merge(cts,gene_ID,by.x="row.names",by.y="ensembl_gene_id")
  Anno.merge.all<-merge(Anno.merge.all,merge.nt,by.x="Row.names",by.y="GeneID")
  names(Anno.merge.all)[1]<-"gene_id" #rename the first column
  Anno.merge.all[Anno.merge.all==""]<-"-" #replace the empty cells to "-"
  Anno.merge.all<-Anno.merge.all[complete.cases(Anno.merge.all[,grep("pvalue|padj",names(Anno.merge.all))]),] # drop pvalue == NA rows; optional
  hybrid_name<-Anno.merge.all$gene_name # add a column of gene_id/name hybrid
  while ("-" %in% hybrid_name){
    hybrid_name[match("-",hybrid_name)] <-
      Anno.merge.all$gene_id[match("-",hybrid_name)]}
  Anno.merge.all<-add_column(Anno.merge.all, hybrid_name, .after = "gene_name")
  write.csv(Anno.merge.all,file=paste0(sub(".tsv$|.txt$","",filename),"_DEG.csv"),row.names=F)
}

#function for non-host
comparison_nonhost<-function(dds,coldata_vs,filename,cts,merge.nt){
  gene_ID_temp<-list()
  for (i in 1:nrow(coldata_vs)){
    group1<-coldata_vs$group1[i]
    group2<-coldata_vs$group2[i]
    comparison <- results(dds, contrast=c("group", group1, group2))
    comparison_f<-as.data.frame(comparison)
    comparison_f<-comparison_f[order(comparison_f$pvalue),]
    system(paste0("mkdir -p"," ",group1,"_vs_",group2)) #for output file structure
    setwd(paste0(group1,"_vs_",group2))
    write.csv(comparison_f,file=paste0(sub(".tsv$|.txt$","",filename),"_",group1,"_vs_",group2,".csv"))
    setwd("../")
    colnames(comparison_f)<-paste0(colnames(comparison_f),".",group1,"_vs_",group2)
    comparison_f <- tibble::rownames_to_column(comparison_f, "NAME")
    gene_ID_temp[[i]]<-comparison_f
  }
  library(dplyr)
  gene_ID <- purrr::reduce(gene_ID_temp,full_join, by = "NAME") #reduce multiple data.frames in a list to a single data.frame
  names(merge.nt)[1]<-"Name" #rename the first column
  Anno.merge.all<-merge(cts,gene_ID,by.x="row.names", by.y="NAME")
  Anno.merge.all<-merge(Anno.merge.all,merge.nt,by.x="Row.names",by.y="Name")
  names(Anno.merge.all)[1]<-"Name" #rename the first column
  Anno.merge.all[Anno.merge.all==""]<-"-" #replace the empty cells to "-"
  Anno.merge.all<-Anno.merge.all[complete.cases(Anno.merge.all[,grep("pvalue|padj",names(Anno.merge.all))]),] # drop pvalue == NA rows
  write.csv(Anno.merge.all,file=paste0(sub(".tsv$|.txt$","",filename),"_DEG.csv"),row.names=F)
}

# apply the comparison function
if (filename == "host_counts.txt"){
  comparison(dds, coldata_vs, filename,gene_ID,cts,merge.nt)
} else {comparison_nonhost(dds, coldata_vs, filename,cts,merge.nt)}


### Plots ###
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr) #for str_trunc function: restrict the showing character length
library(colorspace)
library(RColorBrewer)
# subset normtrans for visualization; consider all groups; to integrate in the comparison function?? DEG=Anno.merge.all
DEG<-read.csv(paste0(sub(".tsv$|.txt$","",filename),"_DEG.csv"),header = T)
# filter by log2 <0.5 & >0.5 and p-value <0.05 for all groups
flt_groups_all<-data.frame()
for (i in 1:length(grep("pvalue",names(DEG)))){
  p<-DEG[DEG[,grep("pvalue",names(DEG))[i]]<0.05 &
            abs(DEG[,grep("log2FoldChange",names(DEG))[i]])>0.5,]
  flt_groups_all<-rbind(flt_groups_all,p)
}
  # flt_groups<-DEG[DEG[,grep("pvalue",names(DEG))[i]]<0.05 &
  #                   abs(DEG[,grep("log2FoldChange",names(DEG))[i]])>0.5,]
  #
  # flt_groups<-flt_groups[complete.cases(flt_groups),]
  # flt_groups_all <- DEG[(abs(DEG[,"log2FoldChange.ART_Young_vs_Ctl_Young"]) > 0.5 |
  #                  abs(DEG[,"log2FoldChange.ART_Old_vs_ART_Young"]) > 0.5) &
  #                  (DEG[,"pvalue.ART_Young_vs_Ctl_Young"]<0.05 |
  #                     DEG[,"pvalue.ART_Old_vs_ART_Young"]<0.05),]

#subset the interested columnes
subset_normtrans<-flt_groups_all[,grep("normtrans|gene_name|gene_id|Name",names(flt_groups_all))]
if (filename == "host_counts.txt"){
  # pass undefined gene name
  while ("-" %in% subset_normtrans$gene_name){
    subset_normtrans$gene_name[match("-",subset_normtrans$gene_name)] <-
      subset_normtrans$gene_id[match("-",subset_normtrans$gene_name)]
  }
  ##check the duplicated rows
  #n_occur <- data.frame(table(subset_normtrans$gene_name))
  #n_occur[n_occur$Freq > 1,]
  ##drop the duplicated rows
  #flt_groups_all<-unique(flt_groups_all)
  # drop the duplicate rows
  subset_normtrans<-subset_normtrans[!duplicated(subset_normtrans$gene_name),]
  #convert column name (gene_name) to row name
  subset_normtrans2<-subset_normtrans[,-1:-2]
  rownames(subset_normtrans2)<-subset_normtrans[,2]
} else {
  # drop the duplicate rows
  subset_normtrans<-subset_normtrans[!duplicated(subset_normtrans$Name),]
  #convert column name (gene_name) to row name
  subset_normtrans2<-subset_normtrans[,-1]
  rownames(subset_normtrans2)<-subset_normtrans[,1]
}
#chop the sample column names
colnames(subset_normtrans2)<-gsub(".normtrans","",colnames(subset_normtrans2))
#heatmap require matrix
transdata <- as.matrix(subset_normtrans2)
# adjust the rowname length
for (i in 1:length(row.names(transdata))){
  l<-nchar(row.names(transdata)[i])
  if (l > 35){
    row.names(transdata)[i]<-str_trunc(row.names(transdata)[i],35)
  }
}
#read the column annotation information (derive from samplesheet.csv)
Anno_col<-as.data.frame(coldata[,2])
rownames(Anno_col)<-coldata[,1]
colnames(Anno_col)<-"Group"
# annotation color selection
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n<-nrow(unique(Anno_col))
cols = gg_color_hue(n)
names(cols)<-as.character(unique(Anno_col)$Group)
anno_colors<-list(Group=cols)

# draw heat map
library(pheatmap)
hp_thumbnail<-pheatmap(transdata, cluster_cols = T, scale="row",
          annotation_col =Anno_col, annotation_colors=anno_colors,
          show_rownames=F,cex=1)
if (nrow(transdata)>100){
  hp<-pheatmap(transdata, cluster_cols = T, scale="row",
                         fontsize_row=5, annotation_col =Anno_col,
               annotation_colors=anno_colors)
} else {
  hp<-pheatmap(transdata, cluster_cols = T, scale="row",
               annotation_col =Anno_col,
               annotation_colors=anno_colors)
}
# hp<-pheatmap(transdata, cluster_cols = T, scale="row",
#              fontsize_row=5,cex=1, annotation_col =Anno_col, cex=0.9)

# save the plot
ggsave("heatmap_thumbnail.pdf",plot=hp_thumbnail)
if (nrow(transdata)>100){
  ggsave("heatmap.pdf",plot=hp,
                limitsize = F,
                height = 0.1*nrow(transdata),
                width = 0.8*ncol(transdata))
} else {
  ggsave("heatmap.pdf",plot=hp,height = 0.2*nrow(transdata))
}
## draw PCA
if (filename %in%
    c("humann_genefamilies_Abundance_go_translated.tsv",
      "humann_genefamilies_Abundance_kegg_translated.tsv",
      "host_counts.txt")){
  pcadata<-plotPCA(vsd,intgroup=c("group"), returnData=T)
  percentVar<-round(100*attr(pcadata,"percentVar"))

  ggplot(pcadata, aes(PC1, PC2, shape=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    ggtitle("PCA")
  ggsave("PCA.pdf")

  ggplot(pcadata, aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    ggtitle("PCA")
  ggsave("PCA_color.pdf")

  ggplot(pcadata, aes(PC1, PC2, shape=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text_repel(aes(label=pcadata$name),size=3) +
    coord_fixed() +
    theme_bw() +
    ggtitle("PCA")
  ggsave("PCA_label.pdf")

  ggplot(pcadata, aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text_repel(aes(label=pcadata$name),size=3) +
    coord_fixed() +
    theme_bw() +
    ggtitle("PCA")
  ggsave("PCA_label_color.pdf")
}
# draw PCoA, heatmap, diversity tests for microbiome
if (filename %in% c("bracken_species_all",
                    "bracken_phylum_all",
                    "bracken_genus_all")){
  library(phyloseq)
  library(vegan)
  #Imported biom file into phyloseq object
  biomfilename=paste0(dirname(args[1]),"/bracken_species_all.biom")
  data<-import_biom(biomfilename,parseFunction = parse_taxonomy_default)
  colnames(tax_table(data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
  #Estimated and exported alpha diversity
  #The measures of diversity that aren't totally reliant on singletons, eg. Shannon/Simpson, are valid to use, and users can ignore the warning in phyloseq when calculating those measures.
  data.alpha<-estimate_richness(data)
  write.csv(data.alpha, file="alpha-diversity.csv")
  theme_set(theme_bw())
  p.alpha <- plot_richness(data,measures = c("Shannon","Simpson"))
  p.alpha
  ggsave("Alpha_diversity_sample.pdf")
  #Estimated and exported beta diversity (Bray-Curtis)
  #to transformation data by vst before beta diversity analysis
  samples4phyloseq<-as.data.frame(coldata$group, row.names = as.character(coldata$sample_name))
  colnames(samples4phyloseq)<-"group"
  sampledata<-sample_data(samples4phyloseq)
  data<-merge_phyloseq(data,sampledata)
  dds1 <- phyloseq_to_deseq2(data, ~ group)
  data1<-data
  vsd1 <- varianceStabilizingTransformation(dds1)
  vsd1.df<-vsd1@assays@data@listData[[1]]
  vsd1.df[vsd1.df < 0.0] <- 0.0 #adjust negative values after vst
  otu_table(data1) <- otu_table(vsd1.df, taxa_are_rows = TRUE)

  braycurtis <- phyloseq::distance(data1, method = "bray")
  BCmat <- as.matrix(braycurtis)
  write.csv(BCmat, file = "braycurtis.csv")
  #Created and exported PCoA values
  braycurtis.pcoa <- ordinate(data1, method = "PCoA", distance = "bray")
  braycurtis.pcoa.export <- as.data.frame(braycurtis.pcoa$vectors, row.names = NULL, optional = FALSE, cut.names = FALSE, col.names = names(braycurtis.pcoa$vectors), fix.empty.names = TRUE, stringsAsFactors = default.stringsAsFactors())
  write.csv(braycurtis.pcoa.export, file="braycurtis-pcoa.csv")

  # add the group information in coldata
  coldata_order<-coldata # to keep the original order as the samplesheet
  coldata_order$order<-1:nrow(coldata_order)
  braycurtis.pcoa.export<-merge(braycurtis.pcoa.export, coldata_order,by.x="row.names",by.y = "sample_name")
  row.names(braycurtis.pcoa.export)<-braycurtis.pcoa.export[,1]
  braycurtis.pcoa.export<-braycurtis.pcoa.export[,-1]

  ## plot PCoA
  b<-braycurtis.pcoa[["values"]][["Relative_eig"]]
  #barplot(b[b>0],names.arg =colnames(braycurtis.pcoa.export))

  library(forcats)
  # reorder the group column following the value of order column
  braycurtis.pcoa.export %>%
    mutate(group = fct_reorder(group, order)) %>%
    ggplot(aes(Axis.1, Axis.2,color=group)) +
    geom_point(size=3) +
    xlab(paste0("PCoA1: ",round(100*b[1]),"% variance")) +
    ylab(paste0("PCoA2: ",round(100*b[2]),"% variance")) +
    geom_text_repel(aes(label=row.names(braycurtis.pcoa.export)),size=3) +
    coord_fixed() +
    theme_bw() +
    ggtitle("Bray-Curtis Distances PCoA")
  
  braycurtis.pcoa.export %>%
    mutate(group = fct_reorder(group, order)) %>%
    ggplot(aes(Axis.1, Axis.2,color=group)) +
    geom_point(size=3) +
    xlab(paste0("PCoA1: ",round(100*b[1]),"% variance")) +
    ylab(paste0("PCoA2: ",round(100*b[2]),"% variance")) +
    geom_text_repel(aes(label=row.names(braycurtis.pcoa.export)),size=3) +
    coord_fixed() +
    theme_bw() +
    ggtitle("Bray-Curtis Distances PCoA")
  ggsave("PCoA_label_color.pdf")
  
  braycurtis.pcoa.export %>%
    mutate(group = fct_reorder(group, order)) %>%
    ggplot(aes(Axis.1, Axis.2,color=group)) +
    geom_point(size=3) +
    xlab(paste0("PCoA1: ",round(100*b[1]),"% variance")) +
    ylab(paste0("PCoA2: ",round(100*b[2]),"% variance")) +
    coord_fixed() +
    theme_bw() +
    ggtitle("Bray-Curtis Distances PCoA")
  ggsave("PCoA_color.pdf")
  
  braycurtis.pcoa.export %>%
    mutate(group = fct_reorder(group, order)) %>%
    ggplot(aes(Axis.1, Axis.2,color=group)) +
    geom_point(size=3) +
    xlab(paste0("PCoA1: ",round(100*b[1]),"% variance")) +
    ylab(paste0("PCoA2: ",round(100*b[2]),"% variance")) +
    geom_text_repel(aes(label=row.names(braycurtis.pcoa.export)),size=3) +
    coord_fixed() +
    theme_bw() +
    ggtitle("Bray-Curtis Distances PCoA")
  ggsave("PCoA_label.pdf")
  
  braycurtis.pcoa.export %>%
    mutate(group = fct_reorder(group, order)) %>%
    ggplot(aes(Axis.1, Axis.2,color=group)) +
    geom_point(size=3) +
    xlab(paste0("PCoA1: ",round(100*b[1]),"% variance")) +
    ylab(paste0("PCoA2: ",round(100*b[2]),"% variance")) +
    coord_fixed() +
    theme_bw() +
    ggtitle("Bray-Curtis Distances PCoA")
  ggsave("PCoA.pdf")

  # anosim test
  data2<-t(as.matrix(data1@otu_table))
  pathotype.anosim <- anosim(data2, braycurtis.pcoa.export$group)
  # plot results
  pdf("ANOSIM.pdf")
  plot(pathotype.anosim,
       main="Analysis of Diversity in Groups",
       xlab="",
       ylab="")
  dev.off()
  # Start writing to an output file
  sink('ANOSIM-analysis-output.txt')
  summary(pathotype.anosim)
  # Stop writing to the file
  sink()

  # plot heatmap (Phyloseq style)
  sampledata<-as.data.frame(coldata[,2])
  row.names(sampledata)<-coldata[,1]
  colnames(sampledata)<-"Groups"
  sam  = sample_data(sampledata)
  sp = otu_table(norm, taxa_are_rows = TRUE)
  physeq<-phyloseq(sp, sam)
  if (nrow(norm)>100){
    ph<-plot_heatmap(physeq,max.label = nrow(norm)) +
      theme (axis.text.y = element_text(size=(2.2*225/nrow(norm))))
  } else {
    ph<-plot_heatmap(physeq)
  }
  ph$labels$y<-"Species"
  print(ph)
  ggsave("Heatmap_all.png")
  ggsave("Heatmap_all.pdf")

  # Add sample data
  tax  = tax_table(data)
  otu  = otu_table(data)
  data_sam <- phyloseq(otu,tax, sam)

  # change color scale for plot_bar
  HowManyPhyla <- length(levels(as.factor(data_sam@tax_table[,2])))
  getPalette = colorRampPalette(brewer.pal(9, "Set2"))
  PhylaPalette = getPalette(HowManyPhyla)

  # phyloseq bar plots
  plot_bar(data, fill="Phylum") +
    scale_fill_manual(values =PhylaPalette)
  ggsave("Bar_phy.pdf")

  plot_bar(data_sam, fill="Phylum",facet_grid = ~Groups) +
    scale_fill_manual(values =PhylaPalette)
  ggsave("Bar_group_phy.pdf")

  # relative abundance bar plot
  data_sam_relabund<-transform_sample_counts(data_sam, function(x) x / sum(x))

  theme_set(theme_grey())

  plot_bar(data_sam_relabund, fill="Phylum") +
    geom_bar(stat="identity", position="stack") +
    labs(x = "", y = "Relative Abundance\n") +
    theme(panel.background = element_blank()) +
    scale_fill_manual(values =PhylaPalette)
  ggsave("Bar_relative_phy.pdf")

  # alpha diversity box plot
  theme_set(theme_bw())
  plot_richness(data_sam,"Groups",measures = c("Shannon","Simpson")) +
    geom_boxplot()
  ggsave("Alpha_diversity.pdf")
  theme_set(theme_grey())

  # alpha diversity comparisons
  data.alpha_g<-merge(data.alpha,coldata,by.x="row.names",by.y=1)
  library(ggpubr)
  my_comparisons<-list()
  for (i in 1:nrow(coldata_vs)){
    my_comparisons[[i]] <- c(coldata_vs$group1[i],coldata_vs$group2[i])
  }
  ggboxplot(data.alpha_g, x = "group", y = "Shannon",
                             color = "group", palette = "jco",
                             add = "jitter") +
    stat_compare_means(comparisons = my_comparisons,
                       method = "t.test") # Add pairwise comparisons p-value
  ggsave("Alpha_diversity_Shannon.pdf")

  ggboxplot(data.alpha_g, x = "group", y = "Simpson",
                             color = "group", palette = "jco",
                             add = "jitter") +
    stat_compare_means(comparisons = my_comparisons,
                       method = "t.test")
  ggsave("Alpha_diversity_Simpson.pdf")
}

## draw volcano and bar plot
for (i in 1:nrow(coldata_vs)){
  group1<-coldata_vs$group1[i]
  group2<-coldata_vs$group2[i]
  pvalue_name<-paste0("pvalue.",group1,"_vs_",group2)
  log2FoldChange_name<-paste0("log2FoldChange.",group1,"_vs_",group2)
  # nt4v<-cbind(DEG[,"gene_name"],DEG[,pvalue_name],
  #          DEG[,log2FoldChange_name])
  if (filename == "host_counts.txt"){
    nt4v<-cbind(DEG["gene_name"],DEG[pvalue_name],
                DEG[log2FoldChange_name])
    colnames(nt4v)<-c("gene_name","pvalue","log2FoldChange")
    # nt4v<-as.data.frame(nt4v) # volcano plot require data.frame
    # nt4v$pvalue<-as.numeric(nt4v$pvalue) #character to numeric
    # nt4v$log2FoldChange<-as.numeric(nt4v$log2FoldChange)
    # pass undefined gene name
    while ("-" %in% nt4v$gene_name){
      nt4v$gene_name[match("-",nt4v$gene_name)] <-
        DEG$gene_id[match("-",DEG$gene_name)]
    }
    # drop duplicated row
    nt4v<-nt4v[!duplicated(nt4v$gene_name),]
    # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
    # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
    # add a column of NAs
    nt4v$diffexpressed <- "NO"
    # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP"
    nt4v$diffexpressed[nt4v$log2FoldChange > 0.5 & nt4v$pvalue < 0.05] <- "UP"
    # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
    nt4v$diffexpressed[nt4v$log2FoldChange< -0.5 & nt4v$pvalue < 0.05] <- "DOWN"
    ## Volcano plotting
    # v<-ggplot(data=nt4v,
    #           aes(x=log2FoldChange,y=-log10(pvalue)),
    #           col=diffexpressed) +
    #   geom_point() +
    #   theme_minimal()
    # # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold
    # v2<-v+geom_vline(xintercept = c(-0.5,0.5), col="red") +
    #   geom_hline(yintercept=-log10(0.05), col="red")
    # # Change point color
    # mycolors <- c("blue", "red", "black")
    # names(mycolors) <- c("DOWN", "UP", "NO")
    # v3 <- v2 + scale_colour_manual(values = mycolors)

    diffexpressed<-nt4v[nt4v$diffexpressed != "NO",]
    # select top 20 DEG of FC
    topFC<-rbind(top_n(diffexpressed,20,log2FoldChange),top_n(diffexpressed,-20,log2FoldChange))
    topFC$label<-topFC$gene_name
    topFC<-topFC[!duplicated(topFC$gene_name),]
    # Create a new column "label" to nt4v, that will contain the name of genes differentially expressed (NA in case they are not)
    nt4v_label<-merge(nt4v,topFC, all=T)
  } else {
    nt4v<-cbind(DEG["Name"],DEG[pvalue_name],
                DEG[log2FoldChange_name])
    colnames(nt4v)<-c("Name","pvalue","log2FoldChange")
    # pass undefined gene name
    while ("-" %in% nt4v$Name){
      nt4v$Name[match("-",nt4v$Name)] <-
        DEG$gene_id[match("-",DEG$Name)]
    }
    # drop duplicated row
    nt4v<-nt4v[!duplicated(nt4v$Name),]
    # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
    # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
    # add a column of NAs
    nt4v$diffexpressed <- "NO"
    # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP"
    nt4v$diffexpressed[nt4v$log2FoldChange > 0.5 & nt4v$pvalue < 0.05] <- "UP"
    # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
    nt4v$diffexpressed[nt4v$log2FoldChange< -0.5 & nt4v$pvalue < 0.05] <- "DOWN"

    diffexpressed<-nt4v[nt4v$diffexpressed != "NO",]
    # select top 20 DEG of FC
    topFC<-rbind(top_n(diffexpressed,20,log2FoldChange),top_n(diffexpressed,-20,log2FoldChange))
    topFC$label<-topFC$Name
    topFC<-topFC[!duplicated(topFC$Name),]
    # Create a new column "label" to nt4v, that will contain the name of genes differentially expressed (NA in case they are not)
    nt4v_label<-merge(nt4v,topFC, all=T)
  }

  ## Volcano plotting
  # plot adding up all layers we have seen so far
  ggplot(data=nt4v_label, aes(x=log2FoldChange, y=-log10(pvalue),
                              col=diffexpressed, label=label)) +
    geom_point(size = 1.5) +
    theme_minimal() +
    geom_text_repel(size=3) +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.5, 0.5), col="red",linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), col="red",linetype="dashed") +
    ggtitle(paste0(group1,"_vs_",group2)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(group1,"_vs_",group2,"/Volcano_",group1,"_vs_",group2,".pdf"))

  ## bar plotting
  bar <- diffexpressed
  if (nrow(bar)!=0){
    # sort log2FoldChange top to down
    if (filename == "host_counts.txt"){
      bar$c <- with(bar,reorder(gene_name,log2FoldChange))
    } else {
      bar$c <- with(bar,reorder(Name,log2FoldChange))
    }
  
    bar_plot<-ggplot(bar,aes(x=log2FoldChange,y=c,fill=pvalue))+
      geom_bar(stat="identity",aes(fill=diffexpressed), width=0.5) +
      #scale_fill_brewer(palette="Blues")+
      scale_fill_manual(name="Expression",
                        labels = c("Down", "Up"),
                        values = c("DOWN"="#00ba38", "UP"="#f8766d")) +
      #scale_fill_gradient2(low=rgb(14,37,56,max=255),high=rgb(54,169,243,max=255), mid=rgb(52,109,157,max=255), midpoint=0.01)+
      #scale_fill_gradient(low=rgb(54,169,243,max=255), high=rgb(14,37,56,max=255))+
      #scale_fill_continuous_sequential(palette = "Blues3", begin=0.4)+
      labs(x="log2FoldChange",y=" ") + #could add ,fill="n=", title = "BMC_HvsBMC_L_DOWN")
      ggtitle(paste0(group1,"_vs_",group2)) +
      #geom_text(aes(x=log2FoldChange+0.6),label="*") +
      theme_bw() +
      scale_y_discrete(breaks=bar$Name,labels=str_trunc(bar$Name,40)) +
      guides(fill = guide_legend(reverse = TRUE)) #reverse the legend to put "Up" in an upper position
    ggsave(paste0(group1,"_vs_",group2,"/Barplot_",group1,"_vs_",group2,".pdf"),height=0.2*length(bar_plot[["data"]][["Name"]]),limitsize = FALSE)
  }
}

## Venn Diagram
if (length(unique(coldata$group))>=2){
  dv.list<-list()
  for (i in unique(coldata$group)){
    vd<-coldata[coldata$group==i,]
    sample_name<-as.character(vd[,1])
    sample_name.norm<-paste0(sample_name,".norm")
    if (filename == "host_counts.txt"){
      vd.df<-DEG[c("gene_name",sample_name.norm)]
    } else {
      vd.df<-DEG[c("Name",sample_name.norm)]
    }
    #Removing rows having all zeros, ignore 1st column with names
    vd.df<-vd.df[rowSums(vd.df[-1])>0,]
    #add to list
    dv.list[[i]]<-vd.df[,1]
  }

  #myCol <- brewer.pal(length(unique(coldata$group)), "Pastel2")

  library(VennDiagram)
  # Don't write log file for VennDiagram
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  #venn.diagram(...)

  venn.diagram(
    x = dv.list,
    category.names = names(dv.list),
    filename = 'venn_diagramm.png',
    output=TRUE,

    # # Output features
     imagetype="png" ,
    # height = 480 ,
    # width = 480 ,
    # resolution = 300,
     compression = "lzw",

    # Circles
    lwd = 1,
    #lty = 'blank',
    col=cols,
    #fill = myCol,
    # keep same colors with heatmap and PCA
    fill = alpha(cols,0.6),

    # Numbers
    #cex = 1,
    fontface = "bold",
    fontfamily = "sans",

    # Set names
    #cat.cex = 1,
    cat.fontface = "bold",
    # cat.default.pos = "outer",
    # cat.pos = c(-27, 27, 135),
    # cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    #rotation = 1
  )
}


# non-host/host reads ratio comparison
if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
  files_h <- list.files(paste0(dirname(args[1]),"/temp"), pattern="^Report_host_.*\\.txt$", full.names=TRUE, recursive=FALSE)
  lh<-c() #generate a empty list
  la<-c()
  for (i in files_h){
    t<-read.table(i,sep="\t", quote= "")
    a<-t[1,1] #get unclassified reads ratio: non-host part
    r<-a/(100-a) #non-host/host reads ratio in kraken2 host report
    fn<-gsub("^Report_host_|\\.txt$","",basename(i)) #grep the sample name
    lh<-c(lh,setNames(r,fn)) #add sample name with value to the list lh; non-host/host reads ratio
    la<-c(la,setNames(a,fn)) #unclassified reads ratio: non-host part
  }
  lah<-cbind(la,lh)
  lah<-merge(lah,coldata,by.x="row.names",by.y="sample_name")
  colnames(lah)<-c("Sample","unclassified reads ratio %","non-host/host reads ratio","Groups")

  # to draw box plot with comparison of unclassified ratio
  my_comparisons<-list()
  for (i in 1:nrow(coldata_vs)){
    my_comparisons[[i]] <- c(coldata_vs$group1[i],coldata_vs$group2[i])
  }
  ggboxplot(lah, x = "Groups", y = "non-host/host reads ratio",
            color = "Groups", palette = "jco",
            add = "jitter") +
    stat_compare_means(comparisons = my_comparisons,
                       method = "t.test") # Add pairwise comparisons p-value
  ggsave("non-host_vs_host_reads_ratio.pdf")

  ggboxplot(lah, x = "Groups", y = "unclassified reads ratio %",
            color = "Groups", palette = "jco",
            add = "jitter") +
    stat_compare_means(comparisons = my_comparisons,
                       method = "t.test") # Add pairwise comparisons p-value
  ggsave("unclassified_reads_ratio.pdf")
}


# preprocess for graphlan
if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
gph.anno<-read.table(paste0(dirname(args[1]),"/graphlan/annot.txt"),
                     header = F,blank.lines.skip=F, fill=TRUE,sep="\t",
                     comment.char = "",
                     col.names=c("V1","V2","V3","V4"))
gph.anno<-gph.anno[!(gph.anno[,1]=="" & gph.anno[,2]=="annotation" & gph.anno[,3]=="") ,] #clean the Unrecognized annotationline

gph.tree<-read.table(paste0(dirname(args[1]),"/graphlan/tree.txt"),
                     header = F, fill=TRUE)

library(stringi)
sp.name<-row.names(cts)
for ( i in 1:nrow(gph.tree)){
  if (length(strsplit(gph.tree[i,1],"\\.")[[1]])>=7){
    gph.sp<-tail(strsplit(gph.tree[i,1],"\\.")[[1]],1) #extract the last word of the species name from graphlan tree
    gph.sp1<-gsub("_","-",gph.sp) # To match the bracken style
    for ( n in 1:length(sp.name)){
      sp.name1<-tail(strsplit(sp.name," ")[[n]],1) #extract the last word of the species name from bracken file
      if (gph.sp1==sp.name1){
        gph.sp2<-head(tail(strsplit(gph.tree[i,1],"\\.")[[1]],2),1) #extract the second last word of the species name from graphlan tree
        sp.name2<-head(strsplit(sp.name," ")[[n]],1) #extract the first word of the species name from bracken file
        if (grepl(sp.name2,gph.sp2)){
          sp.replace<-sp.name[n]
          sp.replace<-gsub(" ","_",sp.replace)
          #gph.tree[i,1]<-sub(paste0(gph.sp,"([^",gph.sp,"]*)$"),paste0(sp.replace,"\\1"), gph.tree[i,1])
          gph.tree[i,1]<-stri_replace_last_fixed(gph.tree[i,1],gph.sp,sp.replace)
        }
      }
    }
  }
}
write.table(gph.tree,paste0(dirname(args[1]),"/graphlan/tree.txt"),
            col.names = F, row.names = F, quote = F)

start_time <- Sys.time()
sp.replace<-gsub(" ","_",sp.name)
for (i in 31:nrow(gph.anno)){
  for ( n in 1:length(sp.name)){
    sp.name1<-tail(strsplit(sp.name," ")[[n]],1)
    sp.name1<-gsub("-","_",sp.name1) # To match the graphlan anno style
    sp.name2<-head(strsplit(sp.name," ")[[n]],1)
    if (length(grep(sp.name1,gph.anno[[1]])) > 1){
      # check if names match in both of species and genus level
      if (gph.anno[i,1]==sp.name1 & TRUE %in% grepl(sp.name2,gph.anno[(i-30):(i-1),1])){
        gph.anno[i,3]<-gsub(" ","_",gph.anno[i,3])
        if (grepl(gph.anno[i,1],gph.anno[i,3])){
          gph.anno[i,3]<-sub(gph.anno[i,1],sp.replace[n],gph.anno[i,3])
        }
        gph.anno[i,1]<-sp.replace[n]
      }
    } else {
      if (gph.anno[i,1]==sp.name1){
        if (grepl(gph.anno[i,1],gph.anno[i,3])){
          gph.anno[i,3]<-sub(gph.anno[i,1],sp.replace[n],gph.anno[i,3])
        }
        gph.anno[i,1]<-sp.replace[n]
      }
    }
  }
}
Sys.time() - start_time

# add rings to graphlan
ring.list<-list() # make a list for ring color of each group
for (i in unique(coldata$group)){
  ring.gp<-coldata[coldata$group==i,]
  sample_name<-as.character(ring.gp[,1])
  ring.df<-cts[sample_name]
  #Removing rows having all zeros
  ring.df<-ring.df[rowSums(ring.df)>0,]
  ring.df<-ring.df[complete.cases(ring.df),] # remove the NA rows
  ring.df<-gsub(" ","_",row.names(ring.df))
  ring.df<-gsub("\\[||\\]","",ring.df)
  #add to list
  ring.list[[i]]<-ring.df
}

group4gph<-list()
white.list<-data.frame()
for (i in 1:nrow(coldata_vs)){
  group1<-coldata_vs$group1[i]
  group2<-coldata_vs$group2[i]
  gp.vs<-paste0(group1,"_vs_",group2)
  gp.vs.p<-paste0("pvalue.",gp.vs)
  gp.vs.fc<-paste0("log2FoldChange.",gp.vs)
  group.p<-DEG[DEG[gp.vs.p]<0.05 & DEG[gp.vs.fc]>0.5,]
  group.gph<-group.p[["Name"]] # extract significant diff. expressed species
  group.gph<-gsub(" ","_",group.gph)
  group.gph<-gsub("\\[||\\]","",group.gph)
  level<-grep(group1,names(anno_colors$Group))
  ring.color<-ring.list[group1]
  ring.color<-ring.color[[1]]
  ring.color<-ring.color[which(!ring.color %in% group.gph)] #exclude the overlap color
  white.list<-append(white.list,data.frame(temp=group.gph)) # add sig. exp. in a white list with a temporary name
  names(white.list)[names(white.list)=="temp"]<-level #rename to level names
  group.gph<-rbind(c("ring_label", level, group1,""),
                   c("ring_label_color", level, anno_colors$Group[group1],""),
                   c("ring_internal_separator_thickness",level,"1",""),
                   c("ring_external_separator_thickness",level,"1",""),
                   cbind(group.gph,
                         rep("ring_color",length(group.gph)),
                         rep(level,length(group.gph)),
                         rep(anno_colors$Group[group1],length(group.gph))),
                   cbind(ring.color,
                         rep("ring_color",length(ring.color)),
                         rep(level,length(ring.color)),
                         rep(anno_colors$Group[group1],length(ring.color))),
                   cbind(ring.color,
                         rep("ring_alpha",length(ring.color)),
                         rep(level,length(ring.color)),
                         rep(0.2,length(ring.color))))

  group.p2<-DEG[DEG[gp.vs.p]<0.05 & DEG[gp.vs.fc]<(-0.5),]
  group.gph2<-group.p2[["Name"]]
  group.gph2<-gsub(" ","_",group.gph2)
  group.gph2<-gsub("\\[||\\]","",group.gph2)
  level2<-grep(group2,names(anno_colors$Group))
  ring.color2<-ring.list[group2]
  ring.color2<-ring.color2[[1]]
  ring.color2<-ring.color2[which(!ring.color2 %in% group.gph2)] #exclude the overlap color
  white.list<-append(white.list,data.frame(temp=group.gph2)) # add sig. exp. in a white list with a temporary name
  names(white.list)[names(white.list)=="temp"]<-level2 #rename to level names
  group.gph2<-rbind(c("ring_label", level2, group2,""),
                    c("ring_label_color", level2, anno_colors$Group[group2],""),
                    c("ring_internal_separator_thickness",level2,"1",""),
                    c("ring_external_separator_thickness",level2,"1",""),
                    cbind(group.gph2,
                          rep("ring_color",length(group.gph2)),
                          rep(level2,length(group.gph2)),
                          rep(anno_colors$Group[group2],length(group.gph2))),
                    cbind(ring.color2,
                          rep("ring_color",length(ring.color2)),
                          rep(level2,length(ring.color2)),
                          rep(anno_colors$Group[group2],length(ring.color2))),
                    cbind(ring.color2,
                          rep("ring_alpha",length(ring.color2)),
                          rep(level2,length(ring.color2)),
                          rep(0.2,length(ring.color2))))
  group4gph[[i]]<-rbind(unname(group.gph),unname(group.gph2))
}
group4gph.table<-do.call(rbind.data.frame,group4gph)
group4gph.table<-unique(group4gph.table)
names(group4gph.table)[1]<-"V1"

# remove "ring_alpha" from the sig. exp. white list
#group4gph.table1<-group4gph.table[!(group4gph.table[,1]==white.list & group4gph.table[,2]=="ring_alpha"),]
#group4gph.table<-group4gph.table[!((group4gph.table[,1] %in% white.list) & group4gph.table[,2]=="ring_alpha"),]
library(reshape2)
white.ls<-melt(white.list,na.rm=T)
for (g in unique(white.ls$L1)){
  group4gph.table<-group4gph.table[!((group4gph.table[,1] %in% white.ls[white.ls$L1==g,1]) &
                                       group4gph.table[,2]=="ring_alpha" &
                                       group4gph.table[,3]==g),]
}

# gatekeeper of species name format between graphlan and bracken
gph.match<-grep(paste(group4gph.table[[1]],collapse="|"), gph.anno[[1]], value=TRUE)
gph.nomatch<-unique(grep(paste(gph.match,collapse="|"),group4gph.table[[1]],value=TRUE, invert=T))
gph.nomatch<-gph.nomatch[which(!gph.nomatch %in%
                                 c("ring_label","ring_label_color",
                                   "ring_color","ring_alpha",
                                   "ring_internal_separator_thickness",
                                   "ring_external_separator_thickness"))]

for ( n in gph.nomatch){
  gph.nomatch1<-tail(strsplit(n,"_")[[1]],1) #extract the last word of the species name from bracken file
  gph.nomatch2<-head(strsplit(n,"_")[[1]],1) #extract the first word of the species name from bracken file
  count4del<-0
  for (m in gph.anno[[1]]){
    # replace only both 1st and last words matched terms
    if (grepl(gph.nomatch1,m) & grepl(gph.nomatch2,m)){
      group4gph.table[[1]]<-gsub(n, m, group4gph.table[[1]])
      count4del<-count4del+1
    }
  }
  if (count4del==0){
    group4gph.table<-group4gph.table[group4gph.table[,1]!=n,]
  }
}

gph.anno1<-rbind(gph.anno,c("total_plotted_degrees","330","",""),
                 c("start_rotation","270","",""),
                 group4gph.table)
write.table(gph.anno1,paste0(dirname(args[1]),"/graphlan/annot.txt"),
            col.names = F, row.names = F, quote = F, sep="\t", na = "")
}


# adjust covariance effect
normtrans_adj <- limma::removeBatchEffect(normtrans, vsd$group)
# make a folder for outputs
dir.create("../halla",recursive = T)
setwd("../halla")
# save normalized & transformed & covariance corrected data for downstream correlation analysis (e.g. halla)
if (filename == "host_counts.txt"){
  DEG_name<-DEG[c("gene_id","hybrid_name")]
  normtrans_adj<-merge(DEG_name,normtrans_adj,by.x="gene_id",by.y="row.names")
  normtrans_adj<-normtrans_adj[!duplicated(normtrans_adj[2]),]
  #convert column name (gene_name) to row name
  normtrans_adj2<-normtrans_adj[,-1:-2]
  rownames(normtrans_adj2)<-normtrans_adj[,2]
  write.table(normtrans_adj2,"Host_gene.txt",sep="\t",quote=F,col.names=NA)
} else if (filename == "humann_genefamilies_Abundance_go_translated.tsv"){
  write.table(normtrans_adj,"Microbiomes_humann_go.txt",sep="\t",quote=F,col.names=NA)
} else if (filename == "humann_genefamilies_Abundance_kegg_translated.tsv"){
  write.table(normtrans_adj,"Microbiomes_humann_kegg.txt",sep="\t",quote=F,col.names=NA)
} else {
  write.table(normtrans_adj,"Microbiomes.txt",sep="\t",quote=F,col.names=NA)
}
setwd("../")



# for (i in nrow(coldata_vs)){
#   group1<-coldata_vs$group1[i]
#   group2<-coldata_vs$group2[i]
#   colnames_pvalue<-paste0("pvalue",".",group1,"_vs_",group2)
#   colnames_FC<-paste0("log2FoldChange",".",group1,"_vs_",group2)
#   flt_all <- DEG[abs(DEG[,colnames_FC]) > 0.5 & DEG[,colnames_pvalue]<0.05,]
#   flt_all <- DEG[abs(DEG[,colnames_FC]) > 0.5 & DEG[,colnames_pvalue]<0.05,]
#   flt_up <- DEG[DEG[,colnames_FC] > 0.5 & DEG[,colnames_pvalue]<0.05,]
#   flt_down <- DEG[DEG[,colnames_FC] < 0.5 & DEG[,colnames_pvalue]<0.05,]
#   subset_normtrans<-flt_all[,grep("normtrans",names(flt_all))]
#   setwd(paste0(group1,"_vs_",group2))
#   write.csv(flt_all,file=paste0(sub(".tsv$|.txt$","",filename),"_DEG_all.csv"),row.names=F)
#   setwd("../")
# }
