##begin to work on figures

source("../../bin/dermalNFData.R")

##figure 4 is the easiest, starting wtiht hat

library(ggplot2)
library(dplyr)
annotes=rna_annotations()


rna_counts=rna_count_matrix(doLogNorm=TRUE)

res=tidyr::gather(as.data.frame(rna_counts),key=Sample,value=LogCounts)
res$Patient=sapply(annotes$patientId[match(res$Sample,annotes$synapseId)],function(x) gsub("CT0*",'',x))
res$Tissue=annotes$tissueId[match(res$Sample,annotes$synapseId)]
res$Gene=rep(rownames(rna_counts),ncol(rna_counts))
#nres=subset(res,Gene=="NF1")
p<-ggplot(res,aes(y=LogCounts,x=Sample))+geom_violin(aes(fill=Patient))+geom_point(data=subset(res,Gene == 'NF1'),shape=4)
pdf('violinPlotOfLogNormDeseqCounts.pdf')
print(p)
dev.off()

rna_counts=rna_count_matrix(doLogNorm=TRUE,minCount=2)

res=tidyr::gather(as.data.frame(rna_counts),key=Sample,value=LogCounts)
res$Patient=sapply(annotes$patientId[match(res$Sample,annotes$synapseId)],function(x) gsub("CT0*",'',x))
res$Tissue=annotes$tissueId[match(res$Sample,annotes$synapseId)]
res$Gene=rep(rownames(rna_counts),ncol(rna_counts))
#nres=subset(res,Gene=="NF1")
p<-ggplot(res,aes(y=LogCounts,x=Sample))+geom_boxplot(aes(fill=Patient))+geom_point(data=subset(res,Gene == 'NF1'),shape=4)
pdf('boxPlotOfLogNormMinCount2DeseqCounts.pdf')
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
pdf('violinPlotOfSizeFactorNormDeseqCounts.pdf')
print(p)
dev.off()

res$NormCounts=log2(res$NormCounts+1)
p<-ggplot(res,aes(y=NormCounts,x=Sample))+geom_violin(aes(fill=Patient))+geom_point(data=subset(res,Gene == 'NF1'),shape=4)
pdf('violinPlotOfLogSizeFactorNormDeseqCounts.pdf')
print(p)
dev.off()

rna_counts=rna_count_matrix(doNorm=TRUE,minCount=2)

res=tidyr::gather(as.data.frame(rna_counts),key=Sample,value=NormCounts)
res$Patient=sapply(annotes$patientId[match(res$Sample,annotes$synapseId)],function(x) gsub("CT0*",'',x))
res$Tissue=annotes$tissueId[match(res$Sample,annotes$synapseId)]
res$Gene=rep(rownames(rna_counts),ncol(rna_counts))
#nres=subset(res,Gene=="NF1")
p<-ggplot(res,aes(y=NormCounts,x=Sample))+geom_boxplot(aes(fill=Patient))+geom_point(data=subset(res,Gene == 'NF1'),shape=4)+scale_y_log10()
pdf('boxPlotOfSizeFactorNormDeseqCounts.pdf')
print(p)
dev.off()

res$NormCounts=log2(res$NormCounts+1)
p<-ggplot(res,aes(y=NormCounts,x=Sample))+geom_boxplot(aes(fill=Patient))+geom_point(data=subset(res,Gene == 'NF1'),shape=4)
pdf('boxPlotOfLogSizeFactorNormMinCount2DeseqCounts.pdf')
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
pdf('violinPlotOfFpkm.pdf')
print(p)
dev.off()

fres=subset(res,FPKM>0.1)
p<-ggplot(res,aes(y=FPKM+1,x=Sample))+geom_boxplot(aes(fill=Patient))+geom_point(data=subset(res,Gene == 'NF1'),shape=4)+scale_y_log10()
pdf('boxPlotOfFpkmMin0.1.pdf')
print(p)
dev.off()

prot_annote=protein_annotations()
allprots=prot_normalized()

pres=tidyr::gather(as.data.frame(allprots[,-1],key=Sample,value=FoldChange))
colnames(pres)<-c("Sample","FoldChange")
pres$Protein=rep(allprots[,1],ncol(allprots)-1)
pres$Patient=sapply(prot_annote$patientId[match(pres$Sample,prot_annote$synapseId)],function(x) gsub("CT0+","",x))
p<-ggplot(pres,aes(y=FoldChange,x=Sample))+geom_boxplot(aes(fill=Patient))+scale_y_log10()
pdf('boxPlotOfProteins.pdf')
print(p)
dev.off()

p<-ggplot(pres,aes(y=FoldChange,x=Sample))+geom_violin(aes(fill=Patient))+scale_y_log10()
pdf('violinPlotOfProteins.pdf')
print(p)
dev.off()

rnaqc='syn5669831'
protqc='syn5669830'
scripturl='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-06-16/figure_4_5_draft.R'
afiles=list.files('.')
countfiles=afiles[grep('Count',afiles)]
fpkmfiles=afiles[grep('Fpkm',afiles)]
protfiles=afiles[grep('Protein',afiles)]

#now first upload the count files with the count matrix used
for(file in countfiles){
  synStore(File(file,parentId=rnaqc),
           used=list(list(entity='syn5051784',wasExecuted=FALSE),
                     list(url=scripturl,wasExecuted=TRUE)))
}

#now the fpkm
for(file in fpkmfiles){
  synStore(File(file,parentId=rnaqc),
           used=list(list(entity='syn5579598',wasExecuted=FALSE),
                     list(url=scripturl,wasExecuted=TRUE)))
}

#and the proteins
for(file in protfiles){
  synStore(File(file,parentId=protqc),
           used=list(list(entity='syn5305003',wasExecuted=FALSE),
                     list(url=scripturl,wasExecuted=TRUE)))

}








