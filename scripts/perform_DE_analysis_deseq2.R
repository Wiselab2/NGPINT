#!/usr/bin/env Rscript
###################################################################################################################################################
# Code for Sanity Check
###################################################################################################################################################
#commandArgs <- function() c("/opt/rit/app/R/3.3.1/lib64/R/bin/exec/R",                       
#                             "--slave",                                                       
#                             "--no-restore",
#                             "--file=perform_DE_analysis_deseq2.R",
#                             "--args",
#                             "/Users/gsfuerst/Box Sync/Iowa State University/PhD Research/Projects/pElicitor/simulations/",
#                            "/Users/gsfuerst/Box Sync/Iowa State University/PhD Research/Projects/pElicitor/Q001S0001B0001_salmon_counts.matrix",
#                            "0.001",  
#                            "S2B11R1S","S2B11R2S","S2B11R3S",
#                           "S2B11R1N","S2B11R2N","S2B11R3N")
###################################################################################################################################################


# Load required libraries
library(DESeq2)

# Define functions
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Setup variables
args <- commandArgs()
out_dir<-args[6]
print(args)
setwd(out_dir)
counts_file<-args[7]
pval_cutoff<-as.double(args[8])
print(args[9:length(args)])
num_of_replicates<-length(args[9:length(args)])/2

raw_count <- read.delim(counts_file, header=FALSE, row.names=1)
countdata <- round(raw_count)
colSums(countdata)/gm_mean(colSums(countdata))
colnames(countdata)<-args[9:length(args)]

# Remove all gene which has 0 value in all sample
all <- apply(countdata, 1, function(x) all(x==0) )
newdata <- countdata[!all,]

countdata <- as.matrix(newdata)
condition <- factor(c(rep("S",num_of_replicates), rep("B",num_of_replicates)))
coldata <- data.frame(row.names=colnames(countdata), condition)


dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds$condition <- relevel(dds$condition, "B")
sizeFactors(dds)<-colSums(countdata)/gm_mean(colSums(countdata))
dds <- DESeq(dds)

res<-results(dds,independentFiltering=FALSE,cooksCutoff = Inf)
#sort(res$padj)
res<-res[!is.na(res$padj), ]
#res<-res[res$padj<=pval_cutoff,]
res<-res[res$log2FoldChange>0,]
res<-data.frame(res)
write.csv(res,file=paste0(counts_file,"_DGE_Indept_filtering_Off_Cooks_Off.csv"))

res<-results(dds)
res<-res[!is.na(res$padj), ]
res<-res[res$log2FoldChange>0,]
res<-data.frame(res)
write.csv(res,file=paste0(counts_file,"_DGE_Indept_filtering_On_Cooks_On.csv"))

norm_counts<-counts(dds, normalized=TRUE)
write.csv(norm_counts,file=paste0(counts_file,"_norm.csv"))


