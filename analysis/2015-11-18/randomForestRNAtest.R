##random forest classifer of patient samples

library(randomForest)
source("../../bin/dermalNFData.R")

rna.annot<-rna_annotations()
rna<-rna_count_matrix(doNorm=TRUE,minCount=5)

pats<-rna.annot$patientId
names(pats)<-rna.annot$synapseId

df<-data.frame(Patients=as.factor(pats[colnames(rna)]),t(rna))

rf<-randomForest(Patients ~. , data=df)
