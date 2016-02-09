##begin to work on figures

source("../../bin/dermalNFData.R")

##figure 4 is the easiest, starting wtiht hat

library(ggplot2)
library(dplyr)

rna_counts=rna_count_matrix(doLogNorm=TRUE)
annotes=rna_annotations()

res=tidyr::gather(as.data.frame(rna_counts),key=Sample,value=LogCounts)
res$Patient=sapply(annotes$patientId[match(res$Sample,annotes$synapseId)],function(x) gsub("CT0*",'',x))
res$Tissue=annotes$tissueId[match(res$Sample,annotes$synapseId)]
res$Gene=rep(rownames(rna_counts),ncol(rna_counts))
#nres=subset(res,Gene=="NF1")
p<-ggplot(res,aes(y=LogCounts,x=Sample))+geom_violin(aes(fill=Patient))+geom_point(data=subset(res,Gene == 'NF1'),shape=4)
png('violinPlotOfLogNormDeseqCounts.png')
print(p)
dev.off()

##now do regular normalization
rna_counts=rna_count_matrix(doNorm=TRUE)

res=tidyr::gather(as.data.frame(rna_counts),key=Sample,value=NormCounts)
res$Patient=sapply(annotes$patientId[match(res$Sample,annotes$synapseId)],function(x) gsub("CT0*",'',x))
res$Tissue=annotes$tissueId[match(res$Sample,annotes$synapseId)]
res$Gene=rep(rownames(rna_counts),ncol(rna_counts))
#nres=subset(res,Gene=="NF1")
p<-ggplot(res,aes(y=NormCounts,x=Sample))+geom_violin(aes(fill=Patient))+geom_point(data=subset(res,Gene == 'NF1'),shape=4)+scale_y_log10()
png('violinPlotOfSizeFactorNormDeseqCounts.png')
print(p)
dev.off()

res$NormCounts=log2(res$NormCounts+1)
p<-ggplot(res,aes(y=NormCounts,x=Sample))+geom_violin(aes(fill=Patient))+geom_point(data=subset(res,Gene == 'NF1'),shape=4)
png('violinPlotOfLogSizeFactorNormDeseqCounts.png')
print(p)
dev.off()

fpkm_annotes=synapseQuery("select sampleID,patientID,tissueID from entity where parentId=='syn5492805'")

rna_fpkm=rna_fpkm_matrix(FALSE)
samples<-sapply(colnames(rna_fpkm),function(x) gsub('^X','',gsub('.','-',x,fixed=T)))
colnames(rna_fpkm)<-samples
res=tidyr::gather(as.data.frame(rna_fpkm),key=Sample,value=FPKM)
#samples<-sapply(res$Sample,function(x) paste(gsub('^X','',gsub('.','_',x,fixed=T)),'_Final.csv',sep=''))
res$Patient=sapply(fpkm_annotes$entity.patientID[match(res$Sample,fpkm_annotes$entity.sampleID)],function(x) gsub("CT0*",'',x))
res$Tissue=fpkm_annotes$entity.tissueID[match(res$Sample,fpkm_annotes$entity.sampleID)]
res$Gene=rep(rownames(rna_fpkm),ncol(rna_fpkm))
#nres=subset(res,Gene=="NF1")
p<-ggplot(res,aes(y=FPKM+1,x=Sample))+geom_violin(aes(fill=Patient))+geom_point(data=subset(res,Gene == 'NF1'),shape=4)+scale_y_log10()
png('violinPlotOfFpkm.png')
print(p)
dev.off()



##figure 3 focuses on SNP data
cnv_annotes=cnv_annotations()

all.cnv=cnv_unprocessed()
lrr<- do.call("rbind", lapply(all.cnv, function(i) as.data.frame(i[,2:3])))
baf<- do.call("rbind", lapply(all.cnv, function(i) as.data.frame(i[,c(2,4)])))

##figure 5 - plot of proteomics values? 


