##evaluate RNA using glm 

source("../../bin/dermalNFData.R")

dat<-rna_count_matrix()

annot=rna_annotations()
pats=annot$patientId
names(pats)<-annot$synapseId
df<-data.frame(t(dat),Patient=pats[colnames(dat)])

gm<-glm(Patient~.,data=df,family=gaussian())