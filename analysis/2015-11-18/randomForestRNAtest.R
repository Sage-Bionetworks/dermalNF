##random forest classifer of patient samples

library(randomForest)
source("../../bin/dermalNFData.R")

rna.annot<-rna_annotations()
rna<-rna_count_matrix(doNorm=TRUE,minCount=5)

pats<-rna.annot$patientId
names(pats)<-rna.annot$synapseId

df<-data.frame(Patients=as.factor(pats[colnames(rna)]),t(rna))

rf<-randomForest(Patients ~. , data=df)

library(pheatmap)
pheatmap(log10(0.001+rna[rownames(rf$importance)[order(rf$importance,decreasing=T)[1:100]],]),annotation_col=data.frame(Patient=pats),cellheight=10,cellwidth=10,file='top100RFpredictingGenes.png')

nzvals=rf$importance[which(rf$importance>0),]
write(names(nzvals)[order(nzvals,decreasing=T)],file='patientPredictiv')