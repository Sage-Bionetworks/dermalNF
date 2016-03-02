##proteomics QC plots - show clustering by sample
source("../../bin/dermalNFData.R")
library(WGCNA)#for clustering


annotes=proteomics_annotations()



plotDendroAndColors(hclust(as.dist(1-cor(seg.by.region,use='p'))),
                    colors=cbind(rainbow(13)[as.numeric(patients[colnames(seg.by.region)])],
                                 c('black','grey')[as.numeric(as.factor(tissue)[colnames(seg.by.region)])]),
                    dendroLabels=paste("Patient",patients[colnames(seg.by.region)]),cex.dendroLabels=0.7,
                    groupLabels=c('Patient','Tissue'),main='Correlation of CNA Segments by region')
