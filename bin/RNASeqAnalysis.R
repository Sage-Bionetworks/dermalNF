library(ggbiplot)

library(synapseClient)
synapseLogin()

source('../../bin/dermalNFData.R')


##first assemble and store RNA-Seq data
tab<-rna_count_matrix(TRUE)

library(DESeq2)
##then do some column normalization
tab<-t(tab)
annots<-rna_annotations()
annots=annots[match(annots$synapseId,colnames(tab)),]
annots=annots[which(!is.na(annots$patientId)),]
pats=data.frame(PatientID=factor(annots$patientId[match(colnames(tab),annots$synapseId)]))
cds<- DESeqDataSetFromMatrix(tab,colData=pats,~PatientID)#now collect proteomics data

sizeFac<-estimateSizeFactors(cds)

normCounts<-tab/sizeFac@colData$sizeFactor


#now do pca
zv<-which(apply(normCounts,1,var)==0)
if(length(zv)>0)
    normCounts<-normCounts[-zv,]

pc=prcomp(t(normCounts),scale.=TRUE)
pdf('RNASeqClustering.pdf')
                                        #pdf('rnaSeqPCA.pdf')
ggbiplot(pc,groups=annots$patientId,var.axes=FALSE,ellipse=TRUE)

nm=normCounts
colnames(nm)<-paste(annots$patientId,annots$tissueID)
plot(hclust(dist(t(nm))))


dev.off()
library(ggplot)
df<-data.frame(t(nm),Patient=annots$patientId)

pdf("RNANormalizedCounts.pdf")
p<-ggplot(df,aes(Patient,KRAS))+geom_boxplot()
print(p)

p<-ggplot(df,aes(Patient,NF1))+geom_boxplot()
print(p)

p<-ggplot(df,aes(Patient,NRAS))+geom_boxplot()
print(p)

p<-ggplot(df,aes(Patient,BRAF))+geom_boxplot()
print(p)
dev.off()


library(pheatmap)
allvars<-apply(nm,1,var,na.rm=T)
alldisp<-apply(nm,1,function(x) var(x)^2/mean(x))
most.d<-nm[order(alldisp,decreasing=T)[1:100],]
pheatmap(log10(most.d+0.00001),cellheight=10,cellwidth=10,annotation_col=data.frame(Patient=sapply(colnames(most.d),function(x) unlist(strsplit(x,split=' '))[1])),clustering_distance_cols='correlation',clustering_method='ward',file='most100Dispersed.png')

most.d<-nm[order(allvars,decreasing=T)[1:100],]
pheatmap(log10(most.d+0.00001),cellheight=10,cellwidth=10,annotation_col=data.frame(Patient=sapply(colnames(most.d),function(x) unlist(strsplit(x,split=' '))[1])),clustering_distance_cols='correlation',clustering_method='ward',file='most100variable.png')
