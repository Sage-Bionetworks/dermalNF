##format expression data for meta network
source("../../bin/dermalNFData.R")

counts <- rna_count_matrix(stored=TRUE,doVoomNorm=FALSE,minCount=2,doLogNorm=TRUE)

write.table(counts,file='varianceStabilized_RnaSeqCounts_over2.tsv',sep='\t',quote=F)

synStore(File('varianceStabilized_RnaSeqCounts_over2.tsv',parentId=''),used=list(list(url='',wasExecuted=TRUE)))