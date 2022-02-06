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
coldata_factor[]<-lapply(coldata_factor, factor)
coldata<-coldata_factor[,1:2]
# update the coldata if metadata is provided
if (length(args) == 5){
  coldata <- read.csv(args[5], header = T, as.is = F)
  coldata[]<-lapply(coldata, factor)
}

if (filename == "host_counts.txt"){
  # make cts(count matrix) has consistent order with samplesheet/metadata
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
  # make cts(count matrix) has consistent order with samplesheet/metadata
  cts<-cts[coldata$sample_name]
  # load the datastructure to DESeq
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design= ~ group + transcriptome_size)
}

# adjust the design if metadata is provided
funNew <- function(x){
    as.formula(paste("~", paste(x, collapse = " + ")))
  }
if (length(args) == 5){
  design(dds)<-funNew(names(coldata)[2:ncol(coldata)])
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
  library("biomaRt")
  host_sp<-read.csv(args[4]) # read a list of supported host species
  # add annotations; try up to 120 times/20 mins if biomaRt server not response
  ensembl <- NULL
  attempt <- 0
  while ( is.null(ensembl) && attempt <= 120){
    try(
      ensembl <- useMart("ensembl",dataset=host_sp[host_sp$Taxon_ID==args[3],2])
    )
    attempt<-attempt+1
    if (attempt > 1){
      print(paste0("Retry to get gene annotations from ensembl: ",attempt," times")) 
    }
    Sys.sleep(10)
  }
  if (is.null(ensembl)){
    stop("Failed to get gene annotations from ensembl server, please try again later")
  }
  names(merge.nt)[1]<-"GeneID"
  genes <- merge.nt$GeneID
  gene_ID <- getBM(filters="ensembl_gene_id",
                   attributes=c("external_gene_name","ensembl_gene_id",
                                "chromosome_name","start_position","end_position",
                                "strand","gene_biotype","description"),
                   values=genes,mart=ensembl)
  names(gene_ID)[names(gene_ID)=="external_gene_name"]<-"gene_name" #rename the first column
  gene_len<-read.table(args[1], row.names=1, sep="\t",header=T, quote="")
  gene_len<-gene_len["Length"]
  colnames(gene_len)<-"gene_length"
  gene_ID<-merge(gene_ID,gene_len,by.x="ensembl_gene_id", by.y="row.names")
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


## MaAsLin2 ##
library("Maaslin2")
system("mkdir -p MaAsLin2_results")
features <- t(cts)
metadata <- coldata
rownames(metadata)<-metadata[,1]
metadata<-subset(metadata,select=-1)
covar<-toString(shQuote(names(coldata)[2:ncol(coldata)]))
for (r in 1:length(unique(coldata_vs[,2]))){
  fit_data <- Maaslin2(features, metadata, paste0('MaAsLin2_results/ref_',unique(coldata_vs[,2])[r]),
                       fixed_effects = cat(covar, "\n"),
                       reference = paste0("group,",unique(coldata_vs[,2])[r]),
                       plot_heatmap = T, plot_scatter = T,
                       cores=4)
}
if (filename == "host_counts.txt"){ #prepare for gene symbol
  system("mkdir -p MaAsLin2_results/gene_symbol")
  DEG<-read.csv(paste0(sub(".tsv$|.txt$","",filename),"_DEG.csv"),header = T)
  df.DEG<-DEG[!duplicated(DEG$hybrid_name),] #remove duplicate gene symbol
  df4features <- data.frame(df.DEG$hybrid_name, df.DEG[names(cts)], row.names=1)
  features.symbol <- t(df4features)
  for (r in 1:length(unique(coldata_vs[,2]))){
    fit_data <- Maaslin2(features.symbol, metadata, paste0('MaAsLin2_results/gene_symbol/ref_',unique(coldata_vs[,2])[r]),
                         fixed_effects = cat(covar, "\n"),
                         reference = paste0("group,",unique(coldata_vs[,2])[r]),
                         plot_heatmap = T, plot_scatter = T,
                         cores=4)
  }
  write("The number of genes may be less due to the duplicated gene symbols being removed.","MaAsLin2_results/gene_symbol/Readme.txt")
}


if (filename == "host_counts.txt"){
  setwd(paste0(dirname(args[1]),"/Host_DEG")) # go back to the Host_DEG folder
}

### Plots ###
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr) #for str_trunc function: restrict the showing character length
library(colorspace)
library(RColorBrewer)
# subset normtrans for visualization; consider all groups; to integrate in the comparison function?? DEG=Anno.merge.all
DEG<-read.csv(paste0(sub(".tsv$|.txt$","",filename),"_DEG.csv"),header = T)

# count gene TPM
if (filename == "host_counts.txt"){
  RPK <- DEG[,as.character(coldata[,1])]/DEG$gene_length
  TPM <- RPK*1000000/sum(RPK)
  TPM<-cbind(DEG[,c("gene_name","hybrid_name")],TPM)
  write.csv(TPM,file=paste0(sub(".tsv$|.txt$","",filename),"_TPM.csv"),row.names=F)
}

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
  #Imported original biom file into phyloseq object
  biomfilename=paste0(dirname(args[1]),"/temp/bracken_species_all0.biom")
  data<-import_biom(biomfilename,parseFunction = parse_taxonomy_default)
  colnames(tax_table(data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
  #Imported host transcriptome size adjusted & normtrans biom file into phyloseq object
  biomfilename.adj=paste0(dirname(args[1]),"/bracken_species_all.biom")
  data.adj<-import_biom(biomfilename.adj,parseFunction = parse_taxonomy_default)
  colnames(tax_table(data.adj)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
  #Estimated and exported alpha diversity
  #The measures of diversity that aren't totally reliant on singletons, eg. Shannon/Simpson, are valid to use, and users can ignore the warning in phyloseq when calculating those measures.
  data.alpha<-estimate_richness(data.adj)
  write.csv(data.alpha, file="alpha-diversity.csv")
  theme_set(theme_bw())
  p.alpha <- plot_richness(data,measures = c("Shannon","Simpson"))
  p.alpha
  ggsave("Alpha_diversity_sample.pdf")
  #Estimated and exported beta diversity (Bray-Curtis)
  #to transformation data by vst before beta diversity analysis
  if (length(args) == 5){
    samples4phyloseq<-as.data.frame(coldata[2:ncol(coldata)], row.names = as.character(coldata$sample_name))
    sampledata<-sample_data(samples4phyloseq)
    data<-merge_phyloseq(data,sampledata)
    dds1 <- phyloseq_to_deseq2(data, funNew(names(coldata)[2:ncol(coldata)]))
  } else {
    samples4phyloseq<-as.data.frame(coldata$group, row.names = as.character(coldata$sample_name))
    colnames(samples4phyloseq)<-"group"
    sampledata<-sample_data(samples4phyloseq)
    data<-merge_phyloseq(data,sampledata)
    dds1 <- phyloseq_to_deseq2(data, ~ group)
  }
  ## ANCOMBC ##
  library("ANCOMBC")
  pseq <- phyloseq::tax_glom(data, taxrank = "Species") # species taxid
  pseq1 <- microbiome::aggregate_taxa(data,"Species") # species name
  
  ancombc_out <- function(pseq,formula){
    ancombc(
      phyloseq = pseq, 
      formula = formula, 
      p_adj_method = "fdr", 
      zero_cut = 0.90, # by default prevalence filter of 10% is applied
      lib_cut = 0, 
      group = "group", 
      struc_zero = TRUE, 
      neg_lb = TRUE, 
      tol = 1e-5, 
      max_iter = 100, 
      conserve = TRUE, 
      alpha = 0.05, 
      global = TRUE
    )
  } 
  
  if (length(args) == 5){
    out <- ancombc_out(pseq,paste(names(coldata)[2:(ncol(coldata)-1)], collapse = " + "))
    out1 <- ancombc_out(pseq1,paste(names(coldata)[2:(ncol(coldata)-1)], collapse = " + "))
  } else {
    out <- ancombc_out(pseq,"group") # taxid
    out1 <- ancombc_out(pseq1,"group") # species name
  }

  res <- out$res # taxid
  res_global = out$res_global
  res1 <- out1$res # species name
  res_global1 = out1$res_global
  
  system("mkdir -p ANCOMBC_results")
  write.csv(res[["diff_abn"]], file="ANCOMBC_results/diff_abundance.csv")
  write.csv(res[["W"]], file="ANCOMBC_results/Test_statistics.csv")
  write.csv(res[["p_val"]], file="ANCOMBC_results/p_value.csv")
  write.csv(res[["q_val"]], file="ANCOMBC_results/q_value.csv")
  write.csv(res_global, file="ANCOMBC_results/Global_test.csv")
  
  system("mkdir -p ANCOMBC_results/with_species_names")
  write.csv(res1[["diff_abn"]], file="ANCOMBC_results/with_species_names/diff_abundance_name.csv")
  write.csv(res1[["W"]], file="ANCOMBC_results/with_species_names/Test_statistics_name.csv")
  write.csv(res1[["p_val"]], file="ANCOMBC_results/with_species_names/p_value_name.csv")
  write.csv(res1[["q_val"]], file="ANCOMBC_results/with_species_names/q_value_name.csv")
  write.csv(res_global1, file="ANCOMBC_results/with_species_names/Global_test_name.csv")
  
  ## heatmap for ANCOMBC results ##
  ANCOMBC_plot <- function(res){
    library(tidyr)
    heat.mw <- res[["W"]] %>% 
      rownames_to_column("taxid") %>% #preserve rownames
      gather(key, value, -taxid)
    heat.md <- res[["diff_abn"]] %>% 
      rownames_to_column("taxid") %>% #preserve rownames
      gather(key, value, -taxid)
    heat.m<-merge(heat.mw,heat.md,by="taxid")
    
    # Arrange the figure
    p <- ggplot(heat.m, aes(x = key.x, y = taxid, fill = value.x))
    p <- p + geom_tile() 
    p <- p + scale_fill_gradientn("value.x", name ="Test statistics",
                                  breaks = seq(from = -2, to = 2, by = 0.5), 
                                  colours = c("darkblue", "blue", "white", "red", "darkred"), 
                                  limits = c(-2,2)) 
    
    # Polish texts
    p <- p + theme(axis.text.x=element_text(angle = 90, hjust=1, face = "italic"),
                   axis.text.y=element_text(size = 8))
    p <- p + xlab("") + ylab("Taxonomy ID")
    
    # Mark the most significant cells with stars
    if (length(unique(heat.m$taxid))>100){
      p <- p + geom_text(data = subset(heat.m, value.y == "TRUE"), 
                         aes(x = key.x, y = taxid, label = "+"), col = "white", size = 2.5)
    } else {
      p <- p + geom_text(data = subset(heat.m, value.y == "TRUE"), 
                         aes(x = key.x, y = taxid, label = "+"), col = "white", size = 3)
    }
    
    if (length(unique(heat.m$taxid))>100){
      ggsave("heatmap_ANCOMBC.pdf",plot=p,
             limitsize = F,
             height = 0.1*length(unique(heat.m$taxid)),
             width = 1.5+0.6*length(unique(heat.m$key.x)))
    } else {
      ggsave("heatmap_ANCOMBC.pdf",plot=p,height = 0.2*length(unique(heat.m$taxid)))
    }
  }
  setwd("ANCOMBC_results")
  ANCOMBC_plot(res)
  
  ANCOMBC_plot_sp <- function(res1){
    library(tidyr)
    heat.mw <- res1[["W"]] %>% 
      rownames_to_column("name") %>% #preserve rownames
      gather(key, value, -name)
    heat.md <- res1[["diff_abn"]] %>% 
      rownames_to_column("name") %>% #preserve rownames
      gather(key, value, -name)
    heat.m<-merge(heat.mw,heat.md,by="name")
    heat.m$name <- sub(".*s__", "", heat.m$name)
    heat.m$name <- stringr::str_trunc(heat.m$name, 31) # truncate sp_name
    # Arrange the figure
    p <- ggplot(heat.m, aes(x = key.x, y = name, fill = value.x))
    p <- p + geom_tile() 
    p <- p + scale_fill_gradientn("value.x", name ="Test statistics",
                                  breaks = seq(from = -2, to = 2, by = 0.5), 
                                  colours = c("darkblue", "blue", "white", "red", "darkred"), 
                                  limits = c(-2,2)) 
    
    # Polish texts
    p <- p + theme(axis.text.x=element_text(angle = 90, hjust=1, face = "italic"),
                   axis.text.y=element_text(size = 8))
    p <- p + xlab("") + ylab("Species name")
    
    # Mark the most significant cells with stars
    if (length(unique(heat.m$name))>100){
      p <- p + geom_text(data = subset(heat.m, value.y == "TRUE"), 
                         aes(x = key.x, y = name, label = "+"), col = "white", size = 2.5)
    } else {
      p <- p + geom_text(data = subset(heat.m, value.y == "TRUE"), 
                         aes(x = key.x, y = name, label = "+"), col = "white", size = 3)
    }
    
    if (length(unique(heat.m$name))>100){
      ggsave("heatmap_ANCOMBC.pdf",plot=p,
             limitsize = F,
             height = 0.1*length(unique(heat.m$name)),
             width = 2.5+0.8*length(unique(heat.m$key.x)))
    } else {
      ggsave("heatmap_ANCOMBC.pdf",plot=p,height = 0.2*length(unique(heat.m$name)))
    }
  }
  
  setwd("with_species_names")
  ANCOMBC_plot_sp(res1)
  setwd("../..")
  
  # continue for beta diversity
  data1<-data
  vsd1 <- varianceStabilizingTransformation(dds1)
  # normalized reads count with host transcriptome size and with avoiding removing variation associated with the other conditions
  if (length(args) == 5){
    mm <- model.matrix(funNew(names(coldata)[2:(ncol(coldata)-1)]), colData(vsd))
  } else {
    mm <- model.matrix(funNew(names(coldata)[2]), colData(vsd))
  }
  vsd1.df <- limma::removeBatchEffect(assay(vsd1), vsd$transcriptome_size, design=mm)
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
  if (names(braycurtis.pcoa.export)[2]=="Axis.2"){
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
  } else {
    write("No Axis.2 was found on PCoA","No_Axis2_on_PCoA.txt")
  }
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
  tax  = tax_table(data.adj)
  otu  = otu_table(data.adj)
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
      scale_y_discrete(breaks=bar[,1],labels=str_trunc(bar[,1],40)) +
      guides(fill = guide_legend(reverse = TRUE)) #reverse the legend to put "Up" in an upper position
    ggsave(paste0(group1,"_vs_",group2,"/Barplot_",group1,"_vs_",group2,".pdf"),height=0.2*length(bar_plot[["data"]][[1]]),limitsize = FALSE)
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
      if (gph.sp1==sp.name1 || gph.sp==sp.name1){
        gph.sp2<-head(tail(strsplit(gph.tree[i,1],"\\.")[[1]],2),1) #extract the second last word of the species name from graphlan tree
        gph.sp6<-strsplit(gph.tree[i,1],"\\.")[[1]][6] #extract the sixth word of the species name from graphlan tree
        sp.name2<-head(strsplit(sp.name," ")[[n]],1) #extract the first word of the species name from bracken file
        if (grepl(sp.name2,gph.sp2) || grepl(sp.name2,gph.sp6)){
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
if (length(args) == 5){
  coldata.n<-coldata
  coldata.n[]<-lapply(coldata.n, as.numeric)
  normtrans_adj <- limma::removeBatchEffect(normtrans, covariates=coldata.n[,2:ncol(coldata.n)])
}
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


### Pathway enrichment for host genes ###
if (filename == "host_counts.txt"){
  setwd("Host_DEG")
  do.db <- host_sp[host_sp$Taxon_ID==args[3],4] # match host taxID with DO.db database
  # check if variable is NULL
  if(is.null(do.db)){
    print("host species is not supported for pathway enrichment yet")
  } else {
    # if package is not installed, install it
    if (!requireNamespace(do.db, quietly = TRUE))
      BiocManager::install(do.db)
    
    library(clusterProfiler)
    library(enrichplot)
    library(ggnewscale)
    library(do.db, character.only = TRUE)
    library(ggplot2)
    
    # function to make plots for GSEA results
    plots4gsea<-function(edb, data, datax, edb0, genelist, group1, group2){
      # save the top10 GSEA results
      if (nrow(data@result)==0){
        write("No enrichment result was found.",paste0(edb,"/No_enrichment_result_found.txt"))
      } else {
        if (nrow(data@result) >= 10){
          gseaplot2(data, geneSetID = 1:10,rel_heights = c(2, 0.5, 1))
          ggsave(paste0(edb,"/Top10_GSEA_GO.pdf"),height = 10.5, width = 8)
        } else {
          gseaplot2(data, geneSetID = 1:nrow(data@result))
          ggsave(paste0(edb,"/GSEA_GO.pdf"),height = (6+0.5*nrow(data@result)))
        }
        # save all plots of the GSEA GO results
        dir.create(paste0(edb,"/GSEA_all"))
        for (g in 1:nrow(data@result)){
          gseaplot2(data, geneSetID = g, title = data$Description[g])
          ggsave(paste0(edb,"/GSEA_all/",gsub("/","_",data$Description[g]),"_GSEA_",edb0,".pdf"))
        }
        # ridgeline plot for expression distribution of GSEA GO result
        ridgeplot(data)
        if (nrow(data@result)<30){
          ggsave(paste0(edb,"/GSEA_",edb0,"_ridgeplots.pdf"),height = 2+0.48*nrow(data@result))
        } else {
          ggsave(paste0(edb,"/GSEA_",edb0,"_ridgeplots.pdf"), height = 12)
        }
        # dot plot
        dotplot(data) + ggtitle("dotplot for GSEA")
        ggsave(paste0(edb,"/GSEA_",edb0,"_dotplot.pdf"),height = 4.8, width = 6)
        #networks
        cnetplot(datax, foldChange=genelist,cex_label_gene = 0.6)
        if (nrow(data@result)<25){
          ggsave(paste0(edb,"/GSEA_",edb0,"_net.pdf"))
        } else {
          ggsave(paste0(edb,"/GSEA_",edb0,"_net.pdf"),scale= 2)
        }
        # tree plot
        datax2 <- pairwise_termsim(datax)
        treeplot(datax2)
        if (nrow(data@result)<30){
          ggsave(paste0(edb,"/GSEA_",edb0,"_tree.pdf"),height=2+0.48*nrow(data@result),limitsize=F)
        } else {
          ggsave(paste0(edb,"/GSEA_",edb0,"_tree.pdf"),height=6)
        }
        # enrichment map
        emapplot(datax2,layout="kk")
        if (nrow(data@result)<25){
          ggsave(paste0(edb,"/GSEA_",edb0,"_map.pdf"))
        } else {
          ggsave(paste0(edb,"/GSEA_",edb0,"_map.pdf"),scale= 2)
        }
        # Heatmap-like functional classification
        heatplot(datax2, foldChange=genelist)
        if (nrow(data@result)<30){
          ggsave(paste0(edb,"/GSEA_",edb0,"_heat.pdf"),height=2+0.24*nrow(data@result),width=0.2*round(max(nchar((data@result$core_enrichment)))/19),limitsize=F)
        } else {
          ggsave(paste0(edb,"/GSEA_",edb0,"_heat.pdf"),height=6.5,width=0.64*round(max(nchar((data@result$core_enrichment[1:30])))/19),limitsize=F)
        }
        # upset plot
        pdf(file= paste0(edb,"/GSEA_",edb0,"_upset.pdf"),height=0.4*nrow(data@result),width=12)
        p.upset<-upsetplot(data, n=nrow(data@result))
        print(p.upset) # save pdf inside a function
        dev.off()
        # bar plot
        bar<-data@result
        bar$c<-with(bar,reorder(Description,-log(pvalue)))
        p.bar<-ggplot(bar,aes(x=-log(pvalue),y=c,fill=enrichmentScore))+
          geom_bar(stat="identity") +
          scale_fill_gradient2(low="blue",high="red") +
          labs(x="-log(p-value)",y="Description", fill="enrichmentScore") +
          ggtitle(paste0(group1,"_vs_",group2)) +
          #geom_text(aes(x=-log(pvalue)+0.6),label=round(bar$enrichmentScore,2),size= 3) + # could add value on bar
          theme_bw() +
          theme(text = element_text(size = 18))+
          scale_y_discrete(breaks=bar$Description,labels=stringr::str_trunc(bar$Description,40))
        ggsave(paste0(edb,"/GSEA_",edb0,"_barplots.pdf"),plot =p.bar,width=10,height=nrow(data@result)/12*5,limitsize=F)
        
        # to save plots for all results
        if (nrow(data@result)>30){
          ridgeplot(data,showCatdatary = nrow(data@result))
          ggsave(paste0(edb,"/GSEA_",edb0,"_ridgeplots_all.pdf"),height = 0.48*nrow(data@result),limitsize=F)
          dotplot(data, showCatdatary=nrow(data@result)) + ggtitle("dotplot for GSEA")
          ggsave(paste0(edb,"/GSEA_",edb0,"_dotplot_all.pdf"),height = 0.48*nrow(data@result), width = 6,limitsize=F)
          cnetplot(datax, foldChange=genelist,cex_label_gene = 0.6, showCatdatary = round(nrow(data@result)/6))
          ggsave(paste0(edb,"/GSEA_",edb0,"_net_all.pdf"),limitsize=F)
          treeplot(datax2,showCatdatary =nrow(data@result),nCluster = round(nrow(data@result)/6))
          ggsave(paste0(edb,"/GSEA_",edb0,"_tree_all.pdf"),height=0.48*nrow(data@result),limitsize=F)
          emapplot(datax2,layout="kk",showCatdatary = nrow(data@result))
          ggsave(paste0(edb,"/GSEA_",edb0,"_map_all.pdf"),scale=nrow(data@result)/12,limitsize=F)
          heatplot(datax2, foldChange=genelist, howCatdatary = nrow(data@result))
          ggsave(paste0(edb,"/GSEA_",edb0,"_heat_all.pdf"),height=1+0.3*nrow(data@result),width=0.2*round(max(nchar((data@result$core_enrichment)))/19),limitsize=F)
        }
      }
    }
    
    # function for KEGG Pathview plots
    pathview.p<-function(kk,ko.db,kegg_gene_list,dir){
      if (nrow(kk@result)==0){
        write("No enrichment result was found on KEGG.","No_KEGG_enrichment_result.txt")
      } else {
        library("pathview")
        for (g in 1:nrow(kk@result)){
          pathview(gene.data  = kegg_gene_list,
                   pathway.id = kk@result[g,1],
                   species    = ko.db,
                   limit      = list(gene=max(abs(kegg_gene_list)), cpd=1))
        }
      }
    }
    
    # function for pathway enrichment by using comparison results between groups
    enrichment <- function(coldata_vs,args1,do.db){
      for (i in 1:nrow(coldata_vs)){
        group1<-coldata_vs$group1[i]
        group2<-coldata_vs$group2[i]
        setwd(paste0(group1,"_vs_",group2)) # go into each comparison folder
        dir.create("GO")
        dir.create("KEGG")
        dir.create("KEGG/Modules")
        # prepare genelist for GSEA
        df.compare<-read.csv(paste0("host_counts_",group1,"_vs_",group2,".csv"),header = T)
        df.compare<-df.compare[df.compare$pvalue<0.05,]
        # feature 1: numeric vector
        genelist<-df.compare$log2FoldChange
        # feature 2: named vector
        names(genelist) = as.character(df.compare$X)
        # feature 3: decreasing order
        genelist = sort(genelist, decreasing = TRUE)
        
        ## GSEA for GO ##
        ego <- gseGO(geneList     = genelist,
                     OrgDb        = do.db,
                     keyType      = 'ENSEMBL',
                     ont          = "ALL",
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)
        
        # save the full table of GSEA GO results
        write.csv(ego@result,"GO/GSEA_GO_results.csv")
        egox <- setReadable(ego, do.db)
        write.csv(egox@result,"GO/GSEA_GO_results_symbol.csv")
        # draw plots for GO GSEA results
        try(plots4gsea("GO",ego,egox,"GO", genelist, group1, group2))
        
        ## GSEA for KEGG ##
        # KEGG pathway gene set enrichment analysis
        ko.db <- host_sp[host_sp$Taxon_ID==args[3],6] # match host taxID with KEGG database
        
        # Convert gene IDs for gseKEGG function
        # Some genes will be lost here because not all IDs will be converted
        ENS2ENT.ids<-bitr(names(genelist), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=do.db)
        # remove duplicate IDS (e.g. one ENSEMBL match with two ENTREZID)
        dedup_ENS2ENT.ids = ENS2ENT.ids[!duplicated(ENS2ENT.ids[c("ENSEMBL")]),]
        # Create a new dataframe df4kegg which has only the genes which were successfully ID converted
        df4kegg = df.compare[df.compare$X %in% dedup_ENS2ENT.ids$ENSEMBL,]
        # Create a new column in df4kegg with the corresponding ENTREZ IDs
        df4kegg$Y = dedup_ENS2ENT.ids$ENTREZID
        # Create a vector of the gene log2FoldChange
        kegg_gene_list <- df4kegg$log2FoldChange
        # Name vector with ENTREZ ids
        names(kegg_gene_list) <- df4kegg$Y
        # omit any NA values 
        kegg_gene_list<-na.omit(kegg_gene_list)
        # sort the list in decreasing order (required for clusterProfiler)
        kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
        
        # KEGG pathway gene set enrichment analysis
        kk.p <- gseKEGG(geneList     = kegg_gene_list,
                        organism     = ko.db,
                        keyType      = 'ncbi-geneid',
                        minGSSize    = 8,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)
        # KEGG (Modules) functional enrichment analysis
        kk.f <- gseMKEGG(geneList = kegg_gene_list,
                         organism = ko.db,
                         keyType  = 'ncbi-geneid',
                         minGSSize = 8,
                         pvalueCutoff = 0.05)
        
        # save the full table of GSEA KEGG results
        write.csv(kk.p@result,"KEGG/GSEA_KEGG_results.csv")
        kkx.p <- setReadable(kk.p, do.db, 'ENTREZID')
        write.csv(kkx.p@result,"KEGG/GSEA_KEGG_results_symbol.csv")
        write.csv(kk.f@result,"KEGG/Modules/GSEA_KEGG_results.csv")
        kkx.f <- setReadable(kk.f, do.db, 'ENTREZID')
        write.csv(kkx.f@result,"KEGG/Modules/GSEA_KEGG_results_symbol.csv")
        
        # draw plots for KEGG GSEA results
        try(plots4gsea("KEGG",kk.p,kkx.p,"KEGG",kegg_gene_list))
        try(plots4gsea("KEGG/Modules",kk.f,kkx.f,"KEGG",kegg_gene_list))
        
        # KEGG pathview plots
        dir.create("KEGG/Pathview") ; setwd("KEGG/Pathview")
        try(pathview.p(kk.p,ko.db,kegg_gene_list))
        setwd(paste0(dirname(args[1]),"/Host_DEG/",group1,"_vs_",group2)) # go back to each comparison folder from Pathview
        dir.create("KEGG/Modules/Pathview") ; setwd("KEGG/Modules/Pathview")
        try(pathview.p(kk.f,ko.db,kegg_gene_list))
        setwd(paste0(dirname(args[1]),"/Host_DEG/",group1,"_vs_",group2)) # go back to each comparison folder from Pathview
        
        setwd("../") # go back to the Host_DEG folder
      }
    }
    
    ## run GSEA enrichment analysis ##
    enrichment(coldata_vs,args[1],do.db)
    
    # function for preparing genelist for biological theme comparison - compareCluster
    BTC <- function(coldata_vs,do.db){
      genelist.ct<-list()
      for (i in 1:nrow(coldata_vs)){
        group1<-coldata_vs$group1[i]
        group2<-coldata_vs$group2[i]
        # prepare genelist for enrichment
        df.btc<-read.csv(paste0(getwd(),"/",group1,"_vs_",group2,"/","host_counts_",group1,"_vs_",group2,".csv"),header = T)
        # filter DEG
        flt_up <- df.btc[df.btc$log2FoldChange > 0.5 & df.btc$pvalue < 0.05,]
        flt_down <- df.btc[df.btc$log2FoldChange < 0.5 & df.btc$pvalue < 0.05,]
        # up-regulated gene name
        genelist.u<-flt_up$X
        # down-regulated gene name
        genelist.d<-flt_down$X
        # Some genes will be lost here because not all IDs will be converted
        ENS2ENT.u<-bitr(genelist.u, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=do.db)
        ENS2ENT.d<-bitr(genelist.d, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=do.db)
        # remove duplicate IDS (e.g. one ENSEMBL match with two ENTREZID)
        dedup_ENS2ENT.u = ENS2ENT.u[!duplicated(ENS2ENT.u[c("ENSEMBL")]),]
        dedup_ENS2ENT.u = dedup_ENS2ENT.u$ENTREZID
        dedup_ENS2ENT.d = ENS2ENT.d[!duplicated(ENS2ENT.d[c("ENSEMBL")]),]
        dedup_ENS2ENT.d = dedup_ENS2ENT.d$ENTREZID
        # make a list
        genelist.c<-list(dedup_ENS2ENT.u,dedup_ENS2ENT.d)
        names(genelist.c) <- c(paste0(group1,"_vs_",group2,"_UP"), paste0(group1,"_vs_",group2,"_DOWN"))
        genelist.ct<-c(genelist.ct,genelist.c)
      }
      return(genelist.ct)
    }
    ## run biological theme comparison ##
    genelist.ct <- BTC(coldata_vs,do.db)
    # GO enrichment comparison
    try(cgo <- compareCluster(genelist.ct, fun = enrichGO, OrgDb=do.db))
    try(cgo <- setReadable(cgo, OrgDb = do.db, keyType="ENTREZID"))
    if (exists("cgo")==T){
      dotplot(cgo, showCategory = nrow(cgo@compareClusterResult)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
      ggsave("biological_theme_comparison_GO.pdf",height = 0.54*nrow(cgo@compareClusterResult), width = 3*length(unique(cgo@compareClusterResult$Cluster)))
      ggsave("biological_theme_comparison_GO_net.pdf",
             plot = cnetplot(cgo,cex_label_gene = 0.6, showCatdatary = round(nrow(cgo@compareClusterResult)/6)),
             limitsize=F)
    }
    # KEGG enrichment comparison
    try(ck <- compareCluster(genelist.ct, fun = enrichKEGG))
    try(ck <- setReadable(ck, OrgDb = do.db, keyType="ENTREZID"))
    if (exists("ck")==T){
      dotplot(ck, showCategory = nrow(ck@compareClusterResult)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
      ggsave("biological_theme_comparison_KEGG.pdf",height = 0.54*nrow(ck@compareClusterResult), width = 3*length(unique(ck@compareClusterResult$Cluster)))
      ggsave("biological_theme_comparison_KEGG_net.pdf",
             plot = cnetplot(ck,cex_label_gene = 0.6, showCatdatary = round(nrow(ck@compareClusterResult)/6)),
             limitsize=F)
    }
  }
}