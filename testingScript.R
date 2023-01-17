# unit testing using "testthat" on toy dataset - Vered Yoskovitz
library(testthat)
source(file="solution.R")

# test createCountsObject()
mydge=createCountsObject("test\\toyDataset.txt", "test\\annotables_grch38.txt", "test\\airway_metadata.txt")
expect_equal(dim(mydge$counts), c(50,8)) # should return nothing


# test CPMNormalize()
mydge_cpm=CPMNormalize(mydge)
cpmTrue=read.csv("test\\toyDatasetCPM.txt",sep='\t')
cpmTrue=as.matrix(cpmTrue)
expect_equal(round(mydge_cpm$counts,3), round(cpmTrue,3)) # should return nothing

# test filterLowlyExpressedGenes()
mydge=filterLowlyExpressedGenes(mydge)
genesPassedLowloExpressionFilter=c("ENSG00000000003", "ENSG00000000419", "ENSG00000000457", "ENSG00000004455",
                                   "ENSG00000000460", "ENSG00000000938", "ENSG00000000971", "ENSG00000001036", 
                                   "ENSG00000001084", "ENSG00000001167", "ENSG00000001460", "ENSG00000001461", 
                                   "ENSG00000001497", "ENSG00000001561", "ENSG00000001617", "ENSG00000001626", 
                                   "ENSG00000001629", "ENSG00000001630", "ENSG00000001631", "ENSG00000002016", 
                                   "ENSG00000002079", "ENSG00000002330", "ENSG00000002549", "ENSG00000002586", 
                                   "ENSG00000002587", "ENSG00000002726", "ENSG00000002745", "ENSG00000002746", 
                                   "ENSG00000002822", "ENSG00000002834", "ENSG00000002919", "ENSG00000002933", 
                                   "ENSG00000003056", "ENSG00000003096", "ENSG00000003137", "ENSG00000003147", 
                                   "ENSG00000003249", "ENSG00000003393", "ENSG00000003400", "ENSG00000003402", 
                                   "ENSG00000003436", "ENSG00000003509", "ENSG00000003756", "ENSG00000003987", 
                                   "ENSG00000003989", "ENSG00000004059", "ENSG00000004139", "ENSG00000004142", 
                                   "ENSG00000004399")
expect_identical(sort(genesPassedLowloExpressionFilter),sort(rownames(mydge$genes))) # should return nothing
expect_identical(sort(genesPassedLowloExpressionFilter),sort(rownames(mydge$counts))) # should return nothing

# test filterLowlyExpressedGenes()
sig=differentialExpressionAnalysis(mydge)
if (dim(sig)[1]==49) {
  print ("Unit test on toy dataset was successful!")
}
