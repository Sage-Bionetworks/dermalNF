##now that we have the mapping we can compare across mrna/protein and mrna/cnv

source("../../bin/crossDataMapping.R")

if(!exists("alldat")){
                                        #just retrieve data once
    alldat<-getPatientMapping()
}

if(!exists('rna_counts')){
                                        #mapped to gene
    rna_counts=rna_count_matrix(doLogNorm=TRUE,minCount=2)
}

if(!exists('proteomics')){
                                        #mapped to gene
    proteomics=prot_normalized()
}


if(!exists('cnv')){
                                        #cnv mapped to gene
    cnv<-cnv_segmented_by_gene()
}


library(ggplot2)

##function that retrieves matched mRNA and protein levels across common patient and genes
matchMrnaProt<-function(){
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
  return(list(rna=rmat,prot=pmat))

}

#compute spearman rank correlation of patients and genes
mrna_prot_comp<-function(doPlot=TRUE){

    dat.mat<-matchMrnaProt()
    rmat=dat.mat$rna
    pmat=dat.mat$prot
    allpats<-rownames(rmat)
    gene.over<-colnames(rmat)
                                        #get patient Cors
    patientCors<-sapply(allpats,function(x) cor(rmat[x,],pmat[x,],method='spearman'))

                                        #get RNA cors
    geneCors<-sapply(gene.over,function(x) cor(rmat[,x],pmat[,x],method='spearman'))

    if(doPlot){
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
    return(list(geneCorrelation=geneCors,patientCorrelation=patientCors))
}


#retrieve matched mRNA and CNV logR ratio for genes and patients that are in common
matchMrnaCnv<-function(){
  ##first get samples for which there are matched RNA and protein
  over=intersect(which(!is.na(alldat$RNASeq)),which(!is.na(alldat$CNV)))
  over.samps<-alldat[over,]
  print(paste('Found',length(over),'samples with RNASeq and CNV data with',length(unique(over.samps$patient)),'patients'))

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

  return(list(rna=rmat,cnv=pmat))
}

#compute spearman correlation
mrna_cnv_comp<-function(doPlot=TRUE){
    dat.mat<-matchMrnaCnv()
    rmat=dat.mat$rna
    pmat=dat.mat$cnv
    allpats<-rownames(rmat)
    gene.over<-colnames(rmat)

  #get patient Cors
  patientCors<-sapply(allpats,function(x) cor(rmat[x,],pmat[x,],method='spearman'))

  #get RNA cors
  geneCors<-sapply(gene.over,function(x) cor(rmat[,x],pmat[,x],method='spearman'))

    if(doPlot){
                                        #names(allcors)<-allpats
        png('mrnaCNVPatientCorrelations.png')
        res<-ggplot(data.frame(SpearmanCorrelation=patientCors)) + geom_histogram(aes(SpearmanCorrelation)) + labs(title=paste('Spearman Correlation of',nrow(over.samps),'samples\nWith mRNA and CNV Data'))
        print(res)
        dev.off()

        png('mrnaCNVGeneCorrelations.png')
        res<-ggplot(data.frame(SpearmanCorrelation=geneCors)) + geom_histogram(aes(SpearmanCorrelation),binwidth=0.05) + labs(title=paste('Spearman Correlation of',length(gene.over),'genes\nWith mRNA and CNV Data'))
        print(res)
        dev.off()
    }


    return(list(geneCorrelation=geneCors,patientCorrelation=patientCors))

}


getMatchedDistributions<-function(mat1,mat2,numiter=100,method='spearman'){
    ##collect two matrices, permute, then compute row and column correlations,
    patientCors<-c()##how well to random patient samples correlate
    geneCors<-c() ##how well to random genes correlate?

    if(numiter>nrow(mat1)*nrow(mat2))
        print("Number of samples is greater than possible combinations, might be over-sampling samples")
    patientCors<-sapply(1:numiter,function(x) cor(mat1[sample(rownames(mat1),1),],mat2[sample(rownames(mat2),1),],method=method))
    if(numiter>ncol(mat1)*ncol(mat2))
        print("Number of samples is greater than possible combinations, might be over-sampling genes")
    geneCors<-sapply(1:numiter,function(x) cor(mat1[,sample(colnames(mat1),1)],mat2[,sample(colnames(mat2),1)],method=method))

    ## for(i in 1:numiter){
    ##                                     #first permute rows
    ##     rowcors=sapply(sample(rownames(mat1),1),function(m1){
    ##         sapply(sample(rownames(mat2),1),function(m2){
    ##             print(paste(m1,'vs',m2))
    ##             cor(mat1[m1,],mat2[m2,],method=method)})})

    ##     colcors<-sapply(sample(colnames(mat1),1),function(m1){
    ##         sapply(sample(colnames(mat2),1),function(m2){
    ##             print(paste(m1,'vs',m2))
    ##             cor(mat1[,m1],mat2[,m2],method=method)})})
    ##     patientCors<-c(patientCors,rowcors)
    ##     geneCors<-c(geneCors,colcors)
    ## }

    return(list(patientCors=patientCors,geneCors=geneCors))

}
