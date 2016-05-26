##format expression data for meta network
source("../../bin/dermalNFData.R")

counts <- rna_count_matrix(stored=TRUE,doVoomNorm=FALSE,minCount=2,doLogNorm=TRUE)

write.table(counts,file='varianceStabilized_RnaSeqCounts_over2.tsv',sep='\t',quote=F)

synStore(File('varianceStabilized_RnaSeqCounts_over2.tsv',parentId='syn4984701'),used=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-03-18/formatRnaMetaNetwork.R',wasExecuted=TRUE)))