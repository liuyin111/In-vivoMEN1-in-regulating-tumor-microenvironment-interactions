rm(list=ls())
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
library(tidyverse)# The tidyverse is an opinionated collection of R packages designed for data science, as_tibble()
library(tibble)  #rownames_to_column()
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(corrplot)
library(RColorBrewer)
library(pheatmap)
library(stringr) #for function str_split()
library(pathview)
library(org.Hs.eg.db)
library(DOSE)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(enrichplot)

#############load data#################################################
setwd('/Users/user/Desktop/')
tcga  = read.table(file="TCGA-PRAD.htseq_counts.tsv", header = T, row.names=1)
View(tcga[1:20,1:100])
tcga = 2^tcga-1
tcga$gene=rownames(tcga)
tcga$gene=str_split(tcga$gene,'[.]',simplify = T)[,1]

View(tcga[1:20,"gene"])

rownames(tcga)=tcga[,"gene"]
which(colnames(tcga)=="gene")  #output colume number 586
tcga = tcga[1:60481,1:585]  ####remove tail, remove the last 7 NA rows and the last gene column######

write.table(tcga,"tcga_COAD_clean.txt",quote=F,)
tcga_clean  = read.table(file="tcga_COAD_clean.txt", header = T, row.names=1)

###############Perform DESeq analysis to normalize MEN1
pick = tcga

#replace '.' in patient ID with '-'
colnames(pick)=chartr(".", "-", colnames(pick))


pick[, 1:585] <- sapply(pick[, 1:585], as.integer)  ##############

group_list = c(paste("MEN", 1:585, sep=""))  
colData <- data.frame(row.names=colnames(pick), group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = pick,
                              colData = colData,
                              design = ~ group_list)
#pre-filter
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]
View(counts(dds))

###########Generate normalized counts###############################################################
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
# %>% What the function does is to pass the left hand side of the operator to the first argument of the right hand side of the operator.
normalized_counts <- counts(dds, normalized=T) %>% data.frame() %>% rownames_to_column(var="gene")
write.table(round(normalized_counts), file="TCGA_normalized_counts_COAD_512", sep="\t", quote=F, col.names=NA)

############ rank based on MEN1, get patient ID for top 10% and bottom 10%
rownames(normalized_counts)=normalized_counts$gene
normalized_counts = normalized_counts[,2:586] #remove fisrt colomn of gene names

normalized_counts_trans= t(normalized_counts)
normalized_counts_trans = as.data.frame(normalized_counts_trans)


pick_rank = normalized_counts_trans[order(normalized_counts_trans$ENSG00000133895),] #ENSEMBL for MEN1: ENSG00000133895
pick_a = pick_rank[-c(51:535),]
pick_id = data.frame(row.names = rownames(pick_a), MEN1= c(rep("MEN1_low",times=50),rep("MEN1_high",times=50)))

tcga<-tcga_clean
tcga = t(tcga)
as.data.frame(tcga)
idx = match(rownames(pick_id),rownames(tcga))
pick_a=tcga[idx,]
pick_a=t(pick_a)
pick_a=as.data.frame(pick_a)

colnames(pick_a) <- c(paste("MEN1_low", 1:50, sep=""),paste("MEN1_high", 1:50, sep=""))

suppressMessages(library(DESeq2))
group_list = c(rep('MEN1_high', times=50), rep('MEN1_low', times=50)) 

###########Generate normalized counts###############################################################
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=T) %>% data.frame() 
write.table(round(normalized_counts), file="normalized_counts_100_samples", sep="\t", quote=F, col.names=NA)

###################DESeq##
dds <-DESeq(dds)
res <- results(dds, contrast=c("group_list","MEN1_low","MEN1_high"))
## Summarize results
summary(res, alpha = 0.01) 
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG_TCGA=as.data.frame(resOrdered)
DEG_TCGA = na.omit(DEG_TCGA)

############ #ensemble ID to gene SYMBOL 
res_tb <- DEG_TCGA %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
res_tb$ENSEMBL = rownames(res_tb) 
df <- bitr(rownames(res_tb), fromType = "ENSEMBL",
           toType = c( "SYMBOL"),
           OrgDb = org.Hs.eg.db)
head(df)
head(res_tb)

diffSig = res_tb[(res_tb$padj < 0.01 & (res_tb$log2FoldChange > 1 | res_tb$log2FoldChange < (-1))),]
diffSig_down = res_tb[(res_tb$padj < 0.01 & (res_tb$log2FoldChange <(-1))),]
diffSig_up = res_tb[(res_tb$padj < 0.01 & (res_tb$log2FoldChange >1)),]

#####output differential genes#######
write.csv(diffSig,"DESeq_diff_sig.csv")
write.csv(diffSig_up,"DESeq_diff_sig_up.csv")
write.csv(diffSig_down,"DESeq_diff_sig_down.csv")
