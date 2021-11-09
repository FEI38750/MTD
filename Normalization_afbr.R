args = commandArgs(trailingOnly=TRUE) # Passing arguments to an R script from bash/shell command lines

# Firstly, to calculate the normalization of reads
library(DESeq2)
library(tibble)
# Read & preprocess the input file cts before run deseq2 analysis
filename<-basename(args[1])
if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
  # read the data file
  #bracken file (eg. bracken_species_all) without quote symbol "\""; mark empty quote as quote=""
  cts<-read.table(args[1],
                  row.names=1, sep="\t",header=T, quote="")
  # Extract columns with count numbers
  cts<-cts[,grepl("*_num", names(cts))]
  # match column name to sample name
  colnames(cts)<-gsub("^Report_|\\.species.bracken_num$|\\.phylum.bracken_num$|.genus.bracken_num$","",colnames(cts))
}

# read the original samplesheet
coldata0 <- read.csv(args[2], header = T, na.strings=c("","NA"))
# extract the contrast information for reference
coldata_vs<- coldata0[c("group1","group2")]
coldata_vs<-coldata_vs[rowSums(is.na(coldata_vs)) == 0,] #remove the NA rows
# read samplesheet as factors (as.is = F) for Deseq2 statistical analysis
coldata_factor <- read.csv(args[2], header = T, as.is = F)
coldata<-coldata_factor[,1:2]

if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
  files_h <- list.files(path=paste0(dirname(args[1]),"/temp"), pattern="^Report_host_.*\\.txt$", full.names=TRUE, recursive=FALSE)
  # to get a list for host
  transcriptome_size<-c() #generate a empty list
  for (i in files_h){
    t<-read.table(i,sep="\t", quote= "")
    total_reads<-t[1,2] + t[2,2] #get total reads abundance of a sample
    fn<-gsub("^Report_host_|\\.txt$","",basename(i)) #grep the sample name
    transcriptome_size<-c(transcriptome_size,setNames(total_reads,fn)) #add sample name with value to the list lh
  }
  transcriptome_size <- log2(transcriptome_size)-mean(log2(transcriptome_size))
  coldata<-merge(coldata,as.data.frame(transcriptome_size), by.x="sample_name",by.y="row.names")
  # make cts(count matrix) has consistent order with samplesheet
  cts<-cts[coldata$sample_name]
  # load the datastructure to DESeq
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design= ~ group + transcriptome_size)
}

# perform the DESeq analysis
dds <- DESeq(dds)

# normalized reads count
norm<-counts(dds,normalized=T)

# read the un-normalized bracken files (tree like)
files<-list.files(path=args[3],full.names=TRUE, recursive=FALSE)

# Apply normalized results to bracken table (tree like); 3rd column
bracken_l <- list()
for (f in files){
  t<-read.table(f,sep="\t", quote= "")
  i <- basename(f)
  t[2]<-0 # clear the 2nd column value
  for (nr in 1:nrow(norm)){
    for (r in 1:nrow(t)){
      if (row.names(norm)[nr]==trimws(t[r,6])){
        t[r,2:3]<-round(norm[nr,i]) # replace with normalized counts
      }
    }
  }
  spaces<-stringr::str_count(t[,6], "\\G ") 
  level_num <- spaces/2 #Determine which level based on number of spaces
  t[7]<-level_num # add level to a new column
  bracken_l[[i]] <-t # save bracken table (tree like) into a list
}

# Function of applying normalized results to bracken table (tree like); 2nd column
Bracken_adj <- function(t){
  tmp_ll<-1
  for (nl in 1:max(t[7])){ # from F level to above
    pre_r_nearest<-0
    pre_r_nearest_up<-0
    for (r in 1:(nrow(t)-1)){
      if (t[r,7] > t[r+1,7]){ # locate the bottom layer/level S of a branch; store in variable r
        cl <- t[r,7]-nl # current level to fill (target level)
        print(paste0("r= ",r,"  cl= ",cl))
        r_nearest <- match(cl,t[(r+1):nrow(t),7]) # get row number that level number matches to the next one
        r_nearest_up <- match(cl,t[(r-1):1,7]) # get row number that level number matches to the upper one
        if (is.na(r_nearest_up)==F){ # to avoid the minus/external levels
          r_nearest_up <- r-r_nearest_up
          if (is.na(r_nearest)==F){
            r_nearest <- r + r_nearest
            if (r_nearest != pre_r_nearest){ # to avoid multiple r share the same target level
              print(paste0("r_nearest=",r_nearest))
              if (r_nearest_up != pre_r_nearest_up){ # to avoid multiple r share the same target level
                print(paste0("r_nearest_up=",r_nearest_up))
                r_section <- t[r_nearest_up:(r_nearest-1),] # get a subset/section from previous S-after level to the next current level (e.g., G/F/O)
                t[r_nearest_up:(r_nearest-1),][t[r_nearest_up:(r_nearest-1),][,7] == cl,2] <- sum(r_section[r_section[,7]==(cl+1),2]) # save to parent level (F)
                print(t[r_nearest_up:(r_nearest-1),][t[r_nearest_up:(r_nearest-1),][,7] == cl,2])
              }
              pre_r_nearest <- r_nearest
              pre_r_nearest_up <- r_nearest_up
            }
            # if (cl %in% t[isUnique(t[,7]),7]){ # for unique section (R)
            #   t[t[,7]==cl,2]<-sum(t[t[,7]==(t[r,7]-(nl-1)),2])}
          }
          if (is.na(r_nearest)==T){ # from current to bottom of the list
            print(paste0("r_nearest is NA, r=",r))
            print(paste0("pre_nearest_up is ",pre_r_nearest_up))
            print(paste0("r_nearest_up is ",r_nearest_up))
            r_section <- t[r_nearest_up:nrow(t),] # get a subset to the end
            t[r_nearest_up:nrow(t),][t[r_nearest_up:nrow(t),][,7] == cl,2] <- sum(r_section[r_section[,7]==(cl+1),2]) # save to parent level (F)
            print(t[r_nearest_up:nrow(t),][t[r_nearest_up:nrow(t),][,7] == cl,2])
          }
        }
        tmp_ll<-r+1 # next branch starting row number
      }
      if (r == (nrow(t)-1)){ # for the bottom layer/level S in the list
        print(paste0("The last S r=",r,"; tmp_ll=" ,tmp_ll))
        r <- r+1 # locate to the bottom layer/level S in the list; store in variable r
        cl <- t[r,7]-nl # current level to fill (target level)
        r_nearest_up <- match(cl,t[(r-1):tmp_ll,7]) # get row number that level number matches to the upper one
        print(paste0("r_nearest_up is ",r_nearest_up))
        if (is.na(r_nearest_up)==F){ # to avoid the minus/external levels
          r_nearest_up <- r-r_nearest_up
          print(paste0("r_nearest_up is ",r_nearest_up))
          r_section <- t[r_nearest_up:nrow(t),] # get a subset to the end
          t[r_nearest_up:nrow(t),][t[r_nearest_up:nrow(t),][,7] == cl,2] <- sum(r_section[r_section[,7]==(cl+1),2]) # save to parent level (F)
          print(paste0("in last branch ",t[r_nearest_up:nrow(t),][t[r_nearest_up:nrow(t),][,7] == cl,2]))
        }
      }
    }
  }
  return(t)
}

# Run the function to apply normalized results to bracken table (tree like); 2nd column
for (i in 1:length(bracken_l)){
  bracken_l[[i]] <- Bracken_adj(bracken_l[[i]])
}

# Function of correcting the percentage according to the normalized results; 1st column
Bracken_adj1 <- function(t){
  for (i in 1:nrow(t)){
    t[i,1]<-round(t[i,2]/t[1,2]*100, digits = 2)
  }
  t<-t[,-7] # remove the temporary/last column
  return(t)
}

# Run the function to apply percentage of normalized results to bracken table (tree like); 1st column
for (i in 1:length(bracken_l)){
  bracken_l[[i]] <- Bracken_adj1(bracken_l[[i]])
}

# Write tables
for (i in 1:length(bracken_l)){
  write.table(bracken_l[[i]],paste0(args[3],"/",names(bracken_l[i])),sep="\t",quote = F,
              row.names = F, col.names = F)
}

print("Adjustment process is done")
