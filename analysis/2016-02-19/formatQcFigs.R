#'
#' Quality control figures
#' 
#' 
#' 
source('../../bin/dermalNFData.R')
require(reshape2)
library(WGCNA)
#' Add in clustering to existing qc
makeCNVClusterPlot<-function(){
  qcfolder='syn5669811'
  segdata<-cnv_segmented()
  seg.by.gene=cnv_segmented_by_gene()
  seg.by.region=cnv_segmented_by_region()
  annote=cnv_annotations()
  tissue=annote$tissueType
  names(tissue)<-annote$synapseId

  patients=sapply(annote$patientId,function(x) gsub('CT0+','',x))
  names(patients)=annote$synapseId
  #plotClusterTreeSamples(t(seg.by.gene),y=as.numeric(as.factor(tissue)),traitLabels='Tissue',dendroLabels=patients)
  #plotDendroAndColors(hclust(as.dist(1-cor(seg.by.gene,use='p'))),colors=rainbow(13)[as.numeric(patients[colnames(seg.by.gene)])])
 # par(mfcol=c(2,1))
  pdf("medianLogRByRegion_correlation.pdf")
  plotDendroAndColors(hclust(as.dist(1-cor(seg.by.region,use='p'))),
                      colors=cbind(rainbow(13)[as.numeric(patients[colnames(seg.by.region)])],
                                   c('black','grey')[as.numeric(as.factor(tissue)[colnames(seg.by.region)])]),
                      dendroLabels=paste("Patient",patients[colnames(seg.by.region)]),cex.dendroLabels=0.7,
                        groupLabels=c('Patient','Tissue'),main='Correlation of CNA Segements by region')
  dev.off()
  png("medianLogRByRegion_correlation.png")
  plotDendroAndColors(hclust(as.dist(1-cor(seg.by.region,use='p'))),
                      colors=cbind(rainbow(13)[as.numeric(patients[colnames(seg.by.region)])],
                                   c('black','grey')[as.numeric(as.factor(tissue)[colnames(seg.by.region)])]),
                      dendroLabels=paste("Patient",patients[colnames(seg.by.region)]),cex.dendroLabels=0.7,
                      groupLabels=c('Patient','Tissue'),main='Correlation of CNA Segements by region')
  dev.off()
  
  synStore(File('medianLogRByRegion_correlation.png',parentId='syn5669811'),
           used=list(list(entity='syn5462067',wasExecuted=FALSE),
                     list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-02-19/formatQcFigs.R',wasExecuted=TRUE)))
           
  synStore(File('medianLogRByRegion_correlation.pdf',parentId='syn5669811'),used=list(list(entity='syn5462067',wasExecuted=FALSE),list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-02-19/formatQcFigs.R',wasExecuted=TRUE)))
  
  
 # leg.text=unique(data.frame(Patient=paste("Patient",patients),Color=rainbow(13)[as.numeric(patients)]))
#  legend("center",legend=leg.text$Patient,col=leg.text$Color)
  }
