source("../../bin/dermalNFData.R")

##first assemble and store RNA-Seq data
tab<-rna_count_files(FALSE)

library(DESeq2)
##then do some column normalization
tab<-t(tab)
annots<-rna_annotations()
pats=data.frame(PatientID=factor(annots$patientId[match(colnames(tab),annots$synapseId)]))
cds<- DESeqDataSetFromMatrix(tab,colData=pats,~PatientID)#now collect proteomics data

sizeFac<-estimateSizeFactors(cds)

normCounts<-tab/sizeFac@colData$sizeFactor

library(ggbiplot)

zv<-which(apply(normCounts,1,var,na.rm=T)==0)
if(length(zv)>0)
  normCounts<-normCounts[-zv,]

pc<-prcomp(t(normCounts),scale.=TRUE,center.=TRUE)
png('normalizedReadCountsClusteredByPatient.png')
ggbiplot(pc,groups=pats[,1],var.axes=F)
dev.off()

#now can we evaluate genes correlated with patient description?
dds<-DESeq(cds)
res<-results(dds)
deg<-rownames(res[which(res$padj<0.05),])

npc<-prcomp(t(normCounts[deg,]),scale.=TRUE,center.=TRUE)
ggbiplot(npc,groups=pats[,1],var.axes=F)


