##proteomics QC plots - show clustering by sample
source("../../bin/dermalNFData.R")
library(WGCNA)#for clustering


#annotes=proteomics_annotations()

#
full.df=prot_unnormalized()
require(ggplot2)

#we can try clustering
require(reshape2)
rat.mat<-acast(full.df,Protein~Sample,value.var='Ratio')
raw.mat<-acast(full.df,Protein~Sample,value.var='RawValue')

zprots=which(apply(raw.mat,1,function(x) all(x==0)))
##now filter by zero proteins

min.df=full.df[-which(full.df$Protein%in%names(zprots)),]

#we can try bbarplots tos how effects of normarlization
p.rat<-ggplot(min.df)+geom_boxplot(aes(x=Tissue,y=Ratio,fill=Experiment))+scale_y_log10()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p.raw<-ggplot(min.df)+geom_boxplot(aes(x=Tissue,y=RawValue,fill=Experiment))+scale_y_log10()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf('proteomicsNormBoxplotsExperiment.pdf')
print(p.rat)
print(p.raw)
dev.off()

p.rat<-ggplot(min.df)+geom_boxplot(aes(x=Tissue,y=Ratio,fill=Patient))+scale_y_log10()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p.raw<-ggplot(min.df)+geom_boxplot(aes(x=Tissue,y=RawValue,fill=Patient))+scale_y_log10()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf('proteomicsNormBoxplotPatient.pdf')
print(p.rat)
print(p.raw)
dev.off()


#

#now do clustering
rat.mat=rat.mat[-zprots,]
raw.mat=raw.mat[-zprots,]


#


pat.exp=t(sapply(unique(full.df$Sample),function(x)
  c(Patient=as.character(full.df$Patient)[match(x,full.df$Sample)],
    Experiment=as.character(full.df$Experiment)[match(x,full.df$Sample)],
    Tissue=as.character(full.df$Tissue)[match(x,full.df$Sample)])))

rownames(pat.exp)<-unique(full.df$Sample)

##first cluster by distance
pdf('distClusters.pdf')
pat.exp=data.frame(pat.exp[colnames(raw.mat),])

plotDendroAndColors(hclust(dist(t(raw.mat)),method='ward'),#as.dist(1-cor(raw.mat,use='p'))),
                    colors=cbind(rainbow(12)[as.numeric(pat.exp$Patient)],
                                 rainbow(14)[as.numeric(pat.exp$Experiment)]),
                    dendroLabels=pat.exp$Tissue,cex.dendroLabels=0.7,
                    groupLabels=c('Patient','Experiment'),main='Correlation of Raw Proteomics Measurements')

pat.exp=data.frame(pat.exp[colnames(rat.mat),])

plotDendroAndColors(hclust(dist(t(rat.mat)),method='ward'),#as.dist(1-cor(rat.mat,use='p'))),
                    colors=cbind(rainbow(12)[as.numeric(pat.exp$Patient)],
                                 rainbow(14)[as.numeric(pat.exp$Experiment)]),
                    dendroLabels=pat.exp$Tissue,cex.dendroLabels=0.7,
                    groupLabels=c('Patient','Experiment'),main='Correlation of Ratio Proteomics Measurements')
dev.off()
##now do log
pdf('logDistClusters.pdf')
pat.exp=data.frame(pat.exp[colnames(raw.mat),])

plotDendroAndColors(hclust(dist(t(log10(1+raw.mat))),method='ward'),#as.dist(1-cor(raw.mat,use='p'))),
                    colors=cbind(rainbow(12)[as.numeric(pat.exp$Patient)],
                                 rainbow(14)[as.numeric(pat.exp$Experiment)]),
                    dendroLabels=pat.exp$Tissue,cex.dendroLabels=0.7,
                    groupLabels=c('Patient','Experiment'),main='Correlation of Raw Proteomics Measurements')

pat.exp=data.frame(pat.exp[colnames(rat.mat),])

plotDendroAndColors(hclust(dist(t(log10(1+rat.mat))),method='ward'),#as.dist(1-cor(rat.mat,use='p'))),
                    colors=cbind(rainbow(12)[as.numeric(pat.exp$Patient)],
                                 rainbow(14)[as.numeric(pat.exp$Experiment)]),
                    dendroLabels=pat.exp$Tissue,cex.dendroLabels=0.7,
                    groupLabels=c('Patient','Experiment'),main='Correlation of Ratio Proteomics Measurements')
#
dev.off()

##then cluster by correlation
pdf('corClusters.pdf')
par(mfrow=c(1,2))
pat.exp=data.frame(pat.exp[colnames(raw.mat),])

plotDendroAndColors(hclust(as.dist(1-cor(raw.mat,use='p')),method='ward'),
                    colors=cbind(rainbow(12)[as.numeric(pat.exp$Patient)],
                                 rainbow(14)[as.numeric(pat.exp$Experiment)]),
                    dendroLabels=pat.exp$Tissue,cex.dendroLabels=0.7,
                    groupLabels=c('Patient','Experiment'),main='Correlation of Raw Proteomics Measurements')

pat.exp=data.frame(pat.exp[colnames(rat.mat),])

plotDendroAndColors(hclust(as.dist(1-cor(rat.mat,use='p')),method='ward'),
                    colors=cbind(rainbow(12)[as.numeric(pat.exp$Patient)],
                                 rainbow(14)[as.numeric(pat.exp$Experiment)]),
                    dendroLabels=pat.exp$Tissue,cex.dendroLabels=0.7,
                    groupLabels=c('Patient','Experiment'),main='Correlation of Ratio Proteomics Measurements')
dev.off()
#now correlation of log
par(mfrow=c(1,2))
pdf("logCorClusters.pdf")
pat.exp=data.frame(pat.exp[colnames(raw.mat),])

plotDendroAndColors(hclust(as.dist(1-cor(log10(1+raw.mat),use='p')),method='ward'),
                    colors=cbind(rainbow(12)[as.numeric(pat.exp$Patient)],
                                 rainbow(14)[as.numeric(pat.exp$Experiment)]),
                    dendroLabels=pat.exp$Tissue,cex.dendroLabels=0.7,
                    groupLabels=c('Patient','Experiment'),main='Correlation of Raw Proteomics Measurements')

pat.exp=data.frame(pat.exp[colnames(rat.mat),])

plotDendroAndColors(hclust(as.dist(1-cor(log10(1+rat.mat),use='p')),method='ward'),
                    colors=cbind(rainbow(12)[as.numeric(pat.exp$Patient)],
                                 rainbow(14)[as.numeric(pat.exp$Experiment)]),
                    dendroLabels=pat.exp$Tissue,cex.dendroLabels=0.7,
                    groupLabels=c('Patient','Experiment'),main='Correlation of Ratio Proteomics Measurements')
dev.off()
