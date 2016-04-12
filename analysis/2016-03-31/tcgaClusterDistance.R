####compare TCGA expression data with dermal data

source('../../bin/TcgaExpressionData.R')
source('../../bin/TcgaMutationalData.R')

##now get dermalNF data

load('combinedExprPC.rda')
library(ggbiplot)
#ggbiplot(pc,var.axes=F)
source('../../bin/dermalNFData.R')
rna.counts<-rna_count_matrix(doLogNorm=T,minCount=2)

tumsByDis$DERMALNF<-colnames(rna.counts)

computeCentroidsForDisease<-function(pc,patlist){
    centroids<-lapply(patlist,function(x){
      pls=intersect(x,rownames(pc$rotation))
      print(paste("Found",length(pls),'samples in PC data'))
      return(colMeans(pc$rotation[pls,]))
    })
    
    names(centroids)<-names(patlist)
    return(centroids)

}

centroids<-computeCentroidsForDisease(pc,tumsByDis)
ggbiplot(prcomp(t(centroids[,-33])),var.axes=F,labels=colnames(centroids)[-33])

dmat<-dist(t(centroids[,-33]))
library(ggplot2)
df<-as.data.frame(as.matrix(dmat)['DERMALNF',])
names(df)<-'Distance'
df$Disease<-factor(rownames(df),levels=rownames(df)[order(df$Distance)])
p<-ggplot(df)+geom_bar(aes(y=Distance,x=Disease),stat='identity')+ theme(axis.text.x=element_text(angle = -90, hjust = 0))

png('distnaceOfDermalToCentroids.png')
print(p)
dev.off()

#sapply(colnames(centroids)[c(1:32,33:35)],function(x) dist(centroids[,x],centroids[,'DERMALNF']))
