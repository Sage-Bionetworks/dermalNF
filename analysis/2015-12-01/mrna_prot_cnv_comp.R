##now that we have the mapping we can compare across mrna/protein and mrna/cnv

source("../../bin/crossDataMapping.R")

#just retrieve data once
alldat<-getPatientMapping()

#mapped to gene
rna_counts=rna_count_matrix(doLogNorm=TRUE,minCount=2)

#mapped to gene
proteomics=prot_normalized()

#cnv mapped to gene
cnv<-cnv_segmented_by_gene()

library(ggplot2)

mrna_prot_comp<-function(){
      
    ##first get samples for which there are matched RNA and protein
  over=intersect(which(!is.na(alldat$RNASeq)),which(!is.na(alldat$Prot)))
  over.samps<-alldat[over,]
  print(paste('Found',length(over),'samples with RNASeq and Proteomics data with',length(unique(over.samps$patient)),'patients'))
  
  gene.over<-intersect(rownames(rna_counts),proteomics$Protein)
  print(paste("Found",length(gene.over),'genes that are in both datasets'))
  
  allcors=c()
  allpats<-c()
  rmat<-NULL
  pmat<-NULL
  for(i in 1:nrow(over.samps)){
    x<-over.samps[i,]
        sampname=paste(x$patient,x$DnaID,x$RnaID,sep='_')
    rvals<-rna_counts[gene.over,x$RNASeq]
    pvals<-proteomics[match(gene.over,proteomics$Protein),x$Prot]
    rmat<-rbind(rmat,rvals)
    pmat<-rbind(pmat,pvals)
    allpats<-c(allpats,sampname)
  
  }
  colnames(rmat)<-colnames(pmat)<-gene.over
  rownames(rmat)<-rownames(pmat)<-allpats
  
  #get patient Cors
  patientCors<-sapply(allpats,function(x) cor(rmat[x,],pmat[x,],method='spearman'))  
  
  #get RNA cors
  geneCors<-sapply(gene.over,function(x) cor(rmat[,x],pmat[,x],method='spearman'))  
  
    png('mrnaProtPatientCorrelations.png')
  res<-ggplot(data.frame(SpearmanCorrelation=patientCors)) + geom_histogram(aes(SpearmanCorrelation),binwidth=0.05) + labs(title=paste('Spearman Correlation of',nrow(over.samps),'samples\nWith mRNA and Proteomics Data'))
  print(res)
  dev.off()
  
  ##now let's figure out which mRNA/prots are correlated? 
  png('mrnaProtGeneCorrelations.png')
  res<-ggplot(data.frame(SpearmanCorrelation=geneCors)) + geom_histogram(aes(SpearmanCorrelation),binwidth=0.05) + labs(title=paste('Spearman Correlation of',length(gene.over),'genes\nWith mRNA and Proteomics Data'))
  print(res)
  dev.off()
}


mrna_cnv_comp<-function(){
  
  ##first get samples for which there are matched RNA and protein
  over=intersect(which(!is.na(alldat$RNASeq)),which(!is.na(alldat$CNV)))
  over.samps<-alldat[over,]
  print(paste('Found',length(over),'samples with RNASeq and Proteomics data with',length(unique(over.samps$patient)),'patients'))
  
  gene.over<-intersect(rownames(rna_counts),rownames(cnv))
  print(paste("Found",length(gene.over),'genes that are in both datasets'))
  
  allcors=c()
  allpats<-c()
  rmat<-NULL
  pmat<-NULL
  for(i in 1:nrow(over.samps)){
    x<-over.samps[i,]
    sampname=paste(x$patient,x$DnaID,x$RnaID,sep='_')
    rvals<-rna_counts[gene.over,x$RNASeq]
    pvals<-cnv[gene.over,x$CNV]
    rmat<-rbind(rmat,rvals)
    pmat<-rbind(pmat,pvals)
    allpats<-c(allpats,sampname)
  }
  colnames(rmat)<-colnames(pmat)<-gene.over
  rownames(rmat)<-rownames(pmat)<-allpats
  #get patient Cors
  patientCors<-sapply(allpats,function(x) cor(rmat[x,],pmat[x,],method='spearman'))  
  
  #get RNA cors
  geneCors<-sapply(gene.over,function(x) cor(rmat[,x],pmat[,x],method='spearman'))  
  
  #names(allcors)<-allpats
  png('mrnaCNVPatientCorrelations.png')
  res<-ggplot(data.frame(SpearmanCorrelation=patientCors)) + geom_histogram(aes(SpearmanCorrelation)) + labs(title=paste('Spearman Correlation of',nrow(over.samps),'samples\nWith mRNA and CNV Data'))
  print(res)
  dev.off()
  
  png('mrnaCNVGeneCorrelations.png')
  res<-ggplot(data.frame(SpearmanCorrelation=geneCors)) + geom_histogram(aes(SpearmanCorrelation),binwidth=0.05) + labs(title=paste('Spearman Correlation of',length(gene.over),'genes\nWith mRNA and CNV Data'))
  print(res)
  dev.off()
  return(list(rna=rmat,cnv=pmat))
  
}