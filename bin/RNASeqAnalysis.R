library(ggbiplot)

library(synapseClient)
synapseLogin()

synq=synapseQuery("select name,id,Patient_ID,Tissue_ID from entity where parentId=='syn4984701'")

synfiles<-sapply(synq$entity.id,synGet)
#now read in all values

allfs<-lapply(synfiles,function(x) read.table(x@filePath,header=T,as.is=T))
names(allfs)<-synq$entity.id

#now get individual genes to create data matrix
hugo.genes<-unique(allfs[[1]][,2])


                                        #now let's get individual counts across patient samples
gene.pat.mat<-sapply(hugo.genes,function(x,allfs){
    res<-sapply(names(allfs),function(y){
        mat<-allfs[[y]]
        sum(mat[which(mat[,2]==x),1])})
    names(res)<-names(allfs)
    res
},allfs)

colnames(gene.pat.mat)<-hugo.genes

##now add in annotations!!!
                                        #now write down matrix to file
write.table(gene.pat.mat,file='GeneBySample.txt',row.names=T,col.names=T)

#now do pca
zv<-which(apply(gene.pat.mat,2,var)==0)
if(length(zv)>0)
    gene.pat.mat<-gene.pat.mat[,-zv]

pc=prcomp(gene.pat.mat,scale.=TRUE)
pdf('RNASeqClustering.pdf')
                                        #pdf('rnaSeqPCA.pdf')
ggbiplot(pc,groups=synq$entity.Patient_ID,var.axes=FALSE)

dev.off()
