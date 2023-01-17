# CytoReason - Exercise for data science candidates, Submitted By Vered Yoskovitz.

library(tidyverse)
library(dplyr)
library(limma)
library(edgeR)
library(ggplot2)
library(limma)

#' @name createCountsObject
#' @description creates an edgeR DGEList object. Genes and samples in the counts matrix that are not annotated are filtered out.
#' Assumptions on the data: 
#' 1. The format is a txt (tab delimited). 
#' 2. The gene names in the count table and in the gene file, and the sample names in the sample file, appear in the first columns.
#' 3. the samples file had the column "type" that assigns samples into groups (for the differential expression analysis)
#' @param countsPath location of counts matrix for the rna-seq data
#' @param genesPath location of gene annotations file for the rna-seq data
#' @param samplesPath location of sample info for the rna-seq data
#' @return edgeR DGEList object, contains all the 3 data files.
#' @example createCountsObject("path\/to\/location\/counts.txt", "path\/to\/location\/gene-annotation.txt","path\/to\/location\/sample-annotation.txt")
createCountsObject <- function(countsPath, genesPath, samplesPath) {
  counts = read.csv(countsPath, sep = '\t')
  genes = read.csv(genesPath, sep = '\t')
  samples = read.csv(samplesPath, sep = '\t')
  names(counts)[1]="geneName"
  names(genes)[1]="geneName"
  genes=genes[!duplicated(genes$geneName), ] # remove duplicated genes
  counts[is.na(counts)] <- 0 # impute missing values with 0
  counts = counts %>% filter(geneName %in% genes[,1]) # keep only genes that are found in the gene-annotation file
  counts = counts[,names(counts) %in% c("geneName", samples[,1])] # keep only samples that are found in the sample-annotation file
  genes=genes %>% filter(geneName %in% counts$geneName) # update the genes matrix to correspond to the updated counts matrix
  
  # make the names as row names for all the 3 dataframes
  rownames(counts)=counts[,1]
  counts=counts[-1]
  rownames(genes)=genes[,1]
  genes=genes[-1]
  rownames(samples)=samples[,1]
  samples=samples[-1]

# create an object  
  dge <-
    DGEList(
      counts = counts,
      genes = genes,
      samples = samples,
      group = samples$type
    )
  return (dge)
}

##################################################
#' @name CPMNormalize
#' @description performs a CPM normalization on the counts matrix
#' @return edgeR DGEList object, where counts attribute is actually a cpm of the input dge object
#' @example CPMNormalize(dge)
CPMNormalize <-function(dgeObj) {
  cpm=sweep(dgeObj$counts, 2, FUN="/", colSums(dgeObj$counts))*10^6 # equivalent to cpm(dgeObj)
  dgeObj$counts=cpm
  return (dgeObj)
}

##################################################
#' @name filterLowlyExpressedGenes
#' @description filters out lowly expressed genes, by the following criteria:
#' keeping genes with CPM>=0.5 in at least 50% of the samples in the group in at least 1 group.
#' In addition, plots the counts of the genes before vs after the filtering.
#' the filtering is performed first on the CPM data, and then applied to the original data.
#' @param dgeObj
#' @return edgeR DGEList object, filtered.
#' @example filterLowlyExpressedGenes(dgeObj)
filterLowlyExpressedGenes <-function(dgeObj) {
  cpmCutoff=0.5
  mingroupsCutoff=1
  genesToKeep=c()
  dgeObj_cpm=CPMNormalize(dgeObj) # apply CPM calculation on the counts data
  originalGeneNumber=dim(dgeObj_cpm$genes)[1]

  for (groupName in unique(dgeObj_cpm$samples$group)) {
    samplesInGroup=dgeObj_cpm$samples %>% filter(group==groupName) %>% rownames()
    minSamplesCutoff=0.5*length(samplesInGroup)
    countsInGroup=dgeObj_cpm$counts[,samplesInGroup]
    numOfSamplesPassedTheGenesCutoff=rowSums((countsInGroup>=cpmCutoff)==TRUE)
    genesToKeepInThisGroup=names(numOfSamplesPassedTheGenesCutoff)[numOfSamplesPassedTheGenesCutoff>=minSamplesCutoff]
    genesToKeep=c(genesToKeep,genesToKeepInThisGroup)
  }
  
  genesToKeep=unique(genesToKeep)
  dgeObj$counts=dgeObj$counts[genesToKeep,]
  dgeObj$genes=dgeObj$genes[genesToKeep,]
  
  geneCount=data.frame(State=c("Before Filtering","After Filtering"),Count=c(originalGeneNumber,dim(dgeObj$genes)[1])) # dataframe for visualization only
  print (ggplot(data=geneCount, aes(x=State,y=Count,fill=State)) + geom_bar(stat="identity")+
    scale_x_discrete(limits=c("Before Filtering","After Filtering")))
  return (dgeObj)
}

##################################################
#' @name differentialExpressionAnalysis
#' @description runs differential Expression Analysis using limma package (this procedure is based on the limma documentation).
#' @return datafrae with all the differential genes (defined by genes with p.adj<=0.05)
#' @example differentialExpressionAnalysis(dge)
differentialExpressionAnalysis <- function(dgeObj){
  if (length(unique(dgeObj$samples$group))<2) {
    print ("DifferentialExpressionAnalysis is not possible")
    return (data.frame())
  }

  mm=model.matrix(~0 + dgeObj$samples$group) # design matrix
  dgeObj <- calcNormFactors(dgeObj) # scale normalization
  v <- voom(dgeObj, design=mm, plot=TRUE)
  fit=lmFit(v,mm)
  fit <- eBayes(fit, trend=TRUE)
  volcanoplot(fit,coef=2,highlight=10,names=dgeObj$genes$SYMBOL)
  allHits=topTable(fit, coef=ncol(mm),number = dim(dgeObj$counts)[1])
  sigGenes=allHits %>% filter(adj.P.Val<=0.05) # filter for significant hits
  return (sigGenes)
}
