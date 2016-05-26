##goal of this script is to visualize mutational data (possibly alongside copy number?)


library(synapseClient)
synapseLogin()

all.muts<-read.table(synGet('syn5611520')@filePath,sep='\t',header=T)
#' Collect germline mutations that are shared across patients
#' 
collectGermlineMuts<-function(minPats=7,maxMuts=20){
    require(reshape2)
  
    gl.muts<-subset(all.muts,Mutation_Status=='Germline')
    ns.gl<-subset(gl.muts,Mutation_Type!='Silent')
    nsi.gl<-subset(ns.gl,Mutation_Type!='Intron')
    mg=names(which(table(nsi.gl$Hugo_Symbol)>minPats))
    mgl=subset(nsi.gl,Hugo_Symbol%in%mg)
    
    mat=acast(mgl,Hugo_Symbol~Sample_ID,value.var='Mutation_Type')
    mat=mat[,order(apply(mat,2,function(x) length(which(x==0))))]
    mat=mat[order(rowSums(mat)),]
    pheatmap(log10(1+mat),cluster_cols=FALSE,cluster_rows=F,cellwidth=10)
    #mat=mat[which(apply(mat,1,max)<maxMuts),]
}

collectSomaticMuts<-function(minPats=1){
  require(dplyr)
  so.muts<-subset(all.muts,Mutation_Status=='Somatic')
  ns.so<-subset(so.muts,Mutation_Type!='Silent')
  nsi.so<-subset(ns.so,Mutation_Type!='Intron')
 
  mg=names(which(table(distinct(nsi.so,Hugo_Symbol,Sample_ID)$Hugo_Symbol)>minPats))
  mgl=subset(nsi.so,Hugo_Symbol%in%mg)
  
  require(reshape2)
  require(pheatmap)
  mat=acast(mgl,Hugo_Symbol~Sample_ID,value.var='Mutation_Type')
  mat=mat[,order(apply(mat,2,function(x) length(which(x==0))))]
  mat=mat[order(rowSums(mat)),]
  pheatmap(log10(1+mat),cluster_cols=F,cluster_rows=F)
  
  pheatmap(log10(1+mat[,which(mat['NF1',]==0)]),cluster_cols=F,cluster_rows=F)
  
  
}