library(ggbiplot)

library(synapseClient)
synapseLogin()



#now do pca
zv<-which(apply(gene.pat.mat,2,var)==0)
if(length(zv)>0)
    gene.pat.mat<-gene.pat.mat[,-zv]

pc=prcomp(gene.pat.mat,scale.=TRUE)
pdf('RNASeqClustering.pdf')
                                        #pdf('rnaSeqPCA.pdf')
ggbiplot(pc,groups=synq$entity.Patient_ID,var.axes=FALSE)

dev.off()
