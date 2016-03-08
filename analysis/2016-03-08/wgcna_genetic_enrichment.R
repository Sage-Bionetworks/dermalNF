##perform enrichment of WGCNA categories
source("../../bin/dermalNFData.R")

require(dplyr)
require(ggplot2)
#signed cluster
#cluster.of.interest<-read.table(synGet('')@filePath,header=T)
#unsigned cluster
cluster.of.interest<-read.table(synGet('syn5714825')@filePath,header=T)
all.mutations<-read.table(synGet('syn5713423')@filePath,header=T,sep='\t')


ann=rna_cufflinks_annotations()
#1. map patient samples to NF1 mutations by position

#2. map patient samples to NF1 mutations by type
#3. map patient samples to cancer germline mutations
#4. map patient samples to cancer somatic mutations (by count)


#4. map patient samples to protein expression